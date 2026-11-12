#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>

#include "rel_geometric_predicates.h"
#include "rel_tetra_mesh_refiner.h"

namespace
{

    //! Number of times a tetrahedron is tried to be refined by splitting the boundary.
    const int maxAttempts = 6;
    //! Number of times the whole mesh is searched for bad tetrahedra.
    const int maxPasses = 4;
    //! Ratio of circumscribed sphere radius to mesh size above which a tetrahedron is too large.
    const double sizeRadiusRatio = 0.75;
    //! A boundary is not split because of a badly shaped tetrahedron, if it would create edges shorter
    //! than this multiple of its shortest edge. A lower value lets edges get shorter and shorter
    //! where parts of the boundary meet in a sharp angle, and the refinement does not end there.
    const double shapeSplitRatio = 1.0;
    //! A boundary is not split because of a too large tetrahedron,
    //! if it would create edges shorter than this multiple of required edge length.
    const double sizeSplitRatio = 0.25;

    const double pi = 3.14159265358979323846;

    inline void subtract(const double *a, const double *b, double *r)
    {
        r[0] = a[0] - b[0];
        r[1] = a[1] - b[1];
        r[2] = a[2] - b[2];
    }

    inline double dot(const double *a, const double *b)
    {
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }

    inline void cross(const double *a, const double *b, double *r)
    {
        r[0] = a[1]*b[2] - a[2]*b[1];
        r[1] = a[2]*b[0] - a[0]*b[2];
        r[2] = a[0]*b[1] - a[1]*b[0];
    }

}

RTetraMeshRefiner::RTetraMeshRefiner(RTetraMesh &mesh, RBoundaryRecovery &boundary, const RTetraMeshOptions &options)
    : mesh(mesh)
    , boundary(boundary)
    , options(options)
    , minLength(0.0)
    , maxSteinerVertices(std::numeric_limits<uint64_t>::max())
    , stopped(false)
{
}

void RTetraMeshRefiner::setMinLength(double minLength)
{
    this->minLength = minLength;
}

void RTetraMeshRefiner::setMaxSteinerVertices(uint64_t maxSteinerVertices)
{
    this->maxSteinerVertices = maxSteinerVertices;
}

bool RTetraMeshRefiner::getStopped() const
{
    return this->stopped;
}

RTetraMeshRefiner::Quality RTetraMeshRefiner::findQuality(int tetraID) const
{
    Quality quality;

    const RTetraMesh::Tetra &rTetra = this->mesh.getTetra(tetraID);
    if (rTetra.region < 0 || this->mesh.isGhost(tetraID))
    {
        return quality;
    }

    const double *p[4] = {this->mesh.getPoint(rTetra.v[0]),
                          this->mesh.getPoint(rTetra.v[1]),
                          this->mesh.getPoint(rTetra.v[2]),
                          this->mesh.getPoint(rTetra.v[3])};

    double u[3], v[3], w[3], vw[3], wu[3], uv[3];
    subtract(p[1],p[0],u);
    subtract(p[2],p[0],v);
    subtract(p[3],p[0],w);
    cross(v,w,vw);
    cross(w,u,wu);
    cross(u,v,uv);

    const double det = dot(u,vw);
    if (!(det > 0.0))
    {
        return quality;
    }
    const double u2 = dot(u,u);
    const double v2 = dot(v,v);
    const double w2 = dot(w,w);

    double r[3];
    for (int i=0;i<3;i++)
    {
        r[i] = (u2 * vw[i] + v2 * wu[i] + w2 * uv[i]) / (2.0 * det);
        quality.center[i] = p[0][i] + r[i];
    }
    quality.radius = std::sqrt(dot(r,r));
    if (!std::isfinite(quality.radius))
    {
        return quality;
    }

    double minEdge2 = std::numeric_limits<double>::max();
    for (int i=0;i<6;i++)
    {
        double edge[3];
        subtract(p[RTetraMesh::edgeVertices[i][1]],p[RTetraMesh::edgeVertices[i][0]],edge);
        minEdge2 = std::min(minEdge2,dot(edge,edge));
    }
    quality.minEdge = std::sqrt(minEdge2);

    // Size.
    const double volume = det / 6.0;
    if (this->options.getVolumeConstraint() > 0.0 && volume > this->options.getVolumeConstraint())
    {
        quality.bad = true;
        quality.tooLarge = true;
        // Edge of a regular tetrahedron of the maximum volume.
        quality.targetEdge = std::cbrt(6.0 * std::sqrt(2.0) * this->options.getVolumeConstraint());
        quality.priority = volume / this->options.getVolumeConstraint();
    }
    if (this->options.getUseSizeFunction())
    {
        double sizeSum = 0.0;
        int nSizes = 0;
        for (int i=0;i<4;i++)
        {
            const double size = this->boundary.getVertexSize(rTetra.v[i]);
            if (size > 0.0)
            {
                sizeSum += size;
                nSizes++;
            }
        }
        if (nSizes > 0)
        {
            const double size = sizeSum / nSizes;
            if (quality.radius > sizeRadiusRatio * size)
            {
                quality.targetEdge = quality.tooLarge ? std::min(quality.targetEdge,size) : size;
                quality.bad = true;
                quality.tooLarge = true;
                quality.priority = std::max(quality.priority,quality.radius / (sizeRadiusRatio * size));
            }
        }
    }
    if (quality.tooLarge)
    {
        // Too large tetrahedra are refined before badly shaped ones.
        quality.priority += 1.0e6;
        return quality;
    }

    // Shape.
    const double ratio = quality.radius / quality.minEdge;
    if (this->options.getRadiusEdgeRatio() > 0.0 && ratio > this->options.getRadiusEdgeRatio())
    {
        quality.bad = true;
        quality.priority = ratio;
        return quality;
    }
    if (this->options.getMinDihedralAngle() > 0.0)
    {
        // Normals of faces pointing inside of the tetrahedron.
        double normals[4][3];
        double lengths[4];
        for (int i=0;i<4;i++)
        {
            const int *pFace = RTetraMesh::faceVertices[i];
            double e1[3], e2[3];
            subtract(p[pFace[1]],p[pFace[0]],e1);
            subtract(p[pFace[2]],p[pFace[0]],e2);
            cross(e1,e2,normals[i]);
            lengths[i] = std::sqrt(dot(normals[i],normals[i]));
        }
        const double limitCos = std::cos(this->options.getMinDihedralAngle() * pi / 180.0);
        double maxCos = -1.0;
        for (int i=0;i<4;i++)
        {
            for (int j=i+1;j<4;j++)
            {
                if (!(lengths[i] > 0.0 && lengths[j] > 0.0))
                {
                    continue;
                }
                const double angleCos = -dot(normals[i],normals[j]) / (lengths[i] * lengths[j]);
                if (angleCos > limitCos)
                {
                    // Faces opposite to vertices i and j share an edge between the other two vertices.
                    // An angle at an edge of the boundary is given by the boundary and can not be made larger.
                    int k = 0;
                    while (k == i || k == j)
                    {
                        k++;
                    }
                    const int l = 6 - i - j - k;
                    if (this->boundary.findSubSegment(rTetra.v[k],rTetra.v[l]) >= 0)
                    {
                        continue;
                    }
                }
                maxCos = std::max(maxCos,angleCos);
            }
        }
        const double minAngle = std::acos(std::min(1.0,std::max(-1.0,maxCos))) * 180.0 / pi;
        if (minAngle < this->options.getMinDihedralAngle())
        {
            quality.bad = true;
            // Lower priority than any tetrahedron with a bad radius/edge ratio.
            quality.priority = 1.0e-3 * (this->options.getMinDihedralAngle() - minAngle);
        }
    }
    return quality;
}

bool RTetraMeshRefiner::refineTetra(const BadTetra &badTetra, const Quality &quality)
{
    if (!(quality.radius > this->minLength))
    {
        return false;
    }

    RBoundaryRecovery::Encroachment encroachment;

    // Walk to the center of circumscribed sphere, which must not lie behind the boundary.
    int currentID = badTetra.tetraID;
    int previousID = RTetraMesh::noTetra;
    int seedID = RTetraMesh::noTetra;
    const uint64_t maxSteps = this->mesh.getNTetras() + 16;
    for (uint64_t step=0;step<maxSteps && encroachment.type==RBoundaryRecovery::EncroachmentNone;step++)
    {
        if (this->mesh.isGhost(currentID))
        {
            return false;
        }
        const RTetraMesh::Tetra &rTetra = this->mesh.getTetra(currentID);
        int nextID = RTetraMesh::noTetra;
        for (int k=0;k<4;k++)
        {
            const int i = (k + int(step)) & 3;
            if (rTetra.n[i] == previousID)
            {
                continue;
            }
            const int *pFace = RTetraMesh::faceVertices[i];
            if (RGeometricPredicates::orient3d(this->mesh.getPoint(rTetra.v[pFace[0]]),
                                               this->mesh.getPoint(rTetra.v[pFace[1]]),
                                               this->mesh.getPoint(rTetra.v[pFace[2]]),
                                               quality.center) < 0.0)
            {
                const int faceID = this->boundary.findSubFace(RTetraMesh::Face{currentID,i});
                if (faceID >= 0)
                {
                    encroachment.type = RBoundaryRecovery::EncroachmentFace;
                    encroachment.id = faceID;
                }
                nextID = rTetra.n[i];
                break;
            }
        }
        if (encroachment.type != RBoundaryRecovery::EncroachmentNone)
        {
            break;
        }
        if (nextID == RTetraMesh::noTetra)
        {
            seedID = currentID;
            break;
        }
        previousID = currentID;
        currentID = nextID;
    }

    if (encroachment.type == RBoundaryRecovery::EncroachmentNone && seedID != RTetraMesh::noTetra)
    {
        // The center may lie right on a subface of the tetrahedron it was found in.
        const RTetraMesh::Tetra &rSeed = this->mesh.getTetra(seedID);
        for (int i=0;i<4;i++)
        {
            const int faceID = this->boundary.findSubFace(RTetraMesh::Face{seedID,i});
            if (faceID < 0)
            {
                continue;
            }
            const int *pFace = RTetraMesh::faceVertices[i];
            if (RGeometricPredicates::orient3d(this->mesh.getPoint(rSeed.v[pFace[0]]),
                                               this->mesh.getPoint(rSeed.v[pFace[1]]),
                                               this->mesh.getPoint(rSeed.v[pFace[2]]),
                                               quality.center) == 0.0)
            {
                encroachment.type = RBoundaryRecovery::EncroachmentFace;
                encroachment.id = faceID;
                break;
            }
        }
    }

    if (encroachment.type == RBoundaryRecovery::EncroachmentNone)
    {
        if (seedID == RTetraMesh::noTetra)
        {
            return false;
        }
        if (this->mesh.findDuplicateVertex(seedID,quality.center) >= 0)
        {
            return false;
        }
        if (!this->boundary.prepareVolumeCavity(quality.center,seedID))
        {
            return false;
        }
        encroachment = this->boundary.findEncroachment(quality.center);
        if (encroachment.type == RBoundaryRecovery::EncroachmentNone)
        {
            this->boundary.insertCavityPoint(quality.center,RBoundaryRecovery::VertexVolume);
            return false;
        }
    }

    // The boundary is in the way, it is split instead. It is not split if the new vertex would be too close
    // to other ones, as shorter edges than those which are there already would be created. Such edges would
    // make new badly shaped tetrahedra, which could never end where parts of the boundary meet in a sharp angle.
    const double minRadius = std::max(2.0 * this->minLength,
                                      quality.tooLarge ? sizeSplitRatio * quality.targetEdge : shapeSplitRatio * quality.minEdge);

    if (encroachment.type == RBoundaryRecovery::EncroachmentSegment)
    {
        return this->boundary.splitSubSegment(encroachment.id,minRadius);
    }
    return this->boundary.splitSubFace(encroachment.id,minRadius);
}

void RTetraMeshRefiner::refine()
{
    this->stopped = false;
    this->boundary.setCollectCreatedTetras(true);

    std::priority_queue<BadTetra> queue;

    auto enqueue = [&](int tetraID, int nAttempts)
    {
        const Quality quality = this->findQuality(tetraID);
        if (!quality.bad)
        {
            return;
        }
        BadTetra badTetra;
        badTetra.priority = quality.priority;
        badTetra.tetraID = tetraID;
        for (int i=0;i<4;i++)
        {
            badTetra.v[i] = this->mesh.getTetra(tetraID).v[i];
        }
        badTetra.nAttempts = nAttempts;
        queue.push(badTetra);
    };

    for (int pass=0;pass<maxPasses && !this->stopped;pass++)
    {
        // Regions are found again, as tetrahedra created at the boundary may have taken a wrong one.
        this->boundary.classifyRegions();

        const int nSlots = this->mesh.getNTetraSlots();
        for (int i=0;i<nSlots;i++)
        {
            if (this->mesh.isAlive(i))
            {
                enqueue(i,0);
            }
        }
        if (queue.empty())
        {
            break;
        }

        const uint64_t nVerticesBefore = this->boundary.getNSteinerVertices();

        while (!queue.empty())
        {
            if (this->boundary.getNSteinerVertices() >= this->maxSteinerVertices)
            {
                this->stopped = true;
                break;
            }

            const BadTetra badTetra = queue.top();
            queue.pop();

            if (!this->mesh.isAlive(badTetra.tetraID))
            {
                continue;
            }
            const RTetraMesh::Tetra &rTetra = this->mesh.getTetra(badTetra.tetraID);
            if (rTetra.v[0] != badTetra.v[0] || rTetra.v[1] != badTetra.v[1] || rTetra.v[2] != badTetra.v[2] || rTetra.v[3] != badTetra.v[3])
            {
                continue;
            }
            const Quality quality = this->findQuality(badTetra.tetraID);
            if (!quality.bad)
            {
                continue;
            }

            const bool tryAgain = this->refineTetra(badTetra,quality);

            // Subsegments and subfaces removed by the new vertex are put back.
            this->boundary.recover();

            std::vector<int> &rCreated = this->boundary.getCreatedTetras();
            for (int tetraID : rCreated)
            {
                if (this->mesh.isAlive(tetraID))
                {
                    enqueue(tetraID,0);
                }
            }
            rCreated.clear();

            if (tryAgain && badTetra.nAttempts + 1 < maxAttempts && this->mesh.isAlive(badTetra.tetraID))
            {
                const RTetraMesh::Tetra &rSame = this->mesh.getTetra(badTetra.tetraID);
                if (rSame.v[0] == badTetra.v[0] && rSame.v[1] == badTetra.v[1] && rSame.v[2] == badTetra.v[2] && rSame.v[3] == badTetra.v[3])
                {
                    BadTetra again = badTetra;
                    again.nAttempts++;
                    queue.push(again);
                }
            }
        }

        while (!queue.empty())
        {
            queue.pop();
        }
        if (this->boundary.getNSteinerVertices() == nVerticesBefore)
        {
            break;
        }
    }

    this->boundary.setCollectCreatedTetras(false);
}
