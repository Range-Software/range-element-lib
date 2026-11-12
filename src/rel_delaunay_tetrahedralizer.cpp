#include <algorithm>
#include <cmath>
#include <cstdint>

#include <rbl_error.h>

#include "rel_delaunay_tetrahedralizer.h"
#include "rel_geometric_predicates.h"

namespace
{

    //! Spread lower 21 bits of a value so that there are two empty bits between each of them.
    inline uint64_t spreadBits(uint64_t value)
    {
        value &= 0x1fffff;
        value = (value | (value << 32)) & 0x1f00000000ffffULL;
        value = (value | (value << 16)) & 0x1f0000ff0000ffULL;
        value = (value | (value << 8)) & 0x100f00f00f00f00fULL;
        value = (value | (value << 4)) & 0x10c30c30c30c30c3ULL;
        value = (value | (value << 2)) & 0x1249249249249249ULL;
        return value;
    }

    inline double findSquaredDistance(const double *p, const double *q)
    {
        return (p[0]-q[0])*(p[0]-q[0]) + (p[1]-q[1])*(p[1]-q[1]) + (p[2]-q[2])*(p[2]-q[2]);
    }

}

void RDelaunayTetrahedralizer::tetrahedralize(RTetraMesh &mesh, std::vector<int> &duplicates)
{
    const int nPoints = mesh.getNPoints();

    duplicates.assign(size_t(nPoints),-1);

    if (nPoints < 4)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,"At least four nodes are needed to generate tetrahedra.");
    }

    // Bounding box.
    double lower[3], upper[3];
    for (int j=0;j<3;j++)
    {
        lower[j] = upper[j] = mesh.getPoint(0)[j];
    }
    for (int i=1;i<nPoints;i++)
    {
        for (int j=0;j<3;j++)
        {
            lower[j] = std::min(lower[j],mesh.getPoint(i)[j]);
            upper[j] = std::max(upper[j],mesh.getPoint(i)[j]);
        }
    }

    // The first tetrahedron is made of points which are far from each other.
    int first[4] = {0,-1,-1,-1};
    {
        double maxValue = 0.0;
        for (int i=0;i<nPoints;i++)
        {
            const double value = findSquaredDistance(mesh.getPoint(i),mesh.getPoint(first[0]));
            if (value > maxValue)
            {
                maxValue = value;
                first[1] = i;
            }
        }
        if (first[1] < 0)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,"All nodes have the same coordinates.");
        }

        const double *p0 = mesh.getPoint(first[0]);
        const double *p1 = mesh.getPoint(first[1]);
        const double u[3] = {p1[0]-p0[0],p1[1]-p0[1],p1[2]-p0[2]};

        maxValue = 0.0;
        for (int i=0;i<nPoints;i++)
        {
            const double *p = mesh.getPoint(i);
            const double w[3] = {p[0]-p0[0],p[1]-p0[1],p[2]-p0[2]};
            const double cx = u[1]*w[2] - u[2]*w[1];
            const double cy = u[2]*w[0] - u[0]*w[2];
            const double cz = u[0]*w[1] - u[1]*w[0];
            const double value = cx*cx + cy*cy + cz*cz;
            if (value > maxValue)
            {
                maxValue = value;
                first[2] = i;
            }
        }
        if (first[2] < 0)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,"All nodes lie on one line.");
        }

        maxValue = 0.0;
        for (int i=0;i<nPoints;i++)
        {
            const double value = std::fabs(RGeometricPredicates::orient3d(p0,p1,mesh.getPoint(first[2]),mesh.getPoint(i)));
            if (value > maxValue)
            {
                maxValue = value;
                first[3] = i;
            }
        }
        if (first[3] < 0)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,"All nodes lie in one plane.");
        }
    }

    mesh.createFirstTetra(first[0],first[1],first[2],first[3]);

    // Points are inserted along a space filling curve, so that each one is close to the previous one.
    struct SortedPoint
    {
        uint64_t key;
        int pointID;
    };
    std::vector<SortedPoint> sortedPoints;
    sortedPoints.reserve(size_t(nPoints));

    double scale[3];
    for (int j=0;j<3;j++)
    {
        const double size = upper[j] - lower[j];
        scale[j] = size > 0.0 ? double(0x1fffff) / size : 0.0;
    }
    for (int i=0;i<nPoints;i++)
    {
        if (i == first[0] || i == first[1] || i == first[2] || i == first[3])
        {
            continue;
        }
        const double *p = mesh.getPoint(i);
        uint64_t key = 0;
        for (int j=0;j<3;j++)
        {
            key |= spreadBits(uint64_t((p[j] - lower[j]) * scale[j])) << j;
        }
        sortedPoints.push_back(SortedPoint{key,i});
    }
    // Points with the same key keep their order, so that a duplicate point follows its original.
    std::stable_sort(sortedPoints.begin(),sortedPoints.end(),[](const SortedPoint &p1, const SortedPoint &p2) { return p1.key < p2.key; });

    for (const SortedPoint &rSortedPoint : sortedPoints)
    {
        int duplicateID = -1;
        if (!mesh.insertVertex(rSortedPoint.pointID,RTetraMesh::noTetra,duplicateID))
        {
            duplicates[size_t(rSortedPoint.pointID)] = duplicateID;
        }
    }
}
