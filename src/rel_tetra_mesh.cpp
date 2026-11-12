#include <algorithm>
#include <cmath>

#include <rbl_error.h>

#include "rel_geometric_predicates.h"
#include "rel_tetra_mesh.h"

const int RTetraMesh::faceVertices[4][3] = {{1,3,2},{0,2,3},{0,3,1},{0,1,2}};

const int RTetraMesh::edgeVertices[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};

namespace
{

    //! Key of an edge which does not depend on order of its vertices.
    inline uint64_t findEdgeKey(int a, int b)
    {
        // Ghost vertex is -1.
        const uint64_t lower = uint64_t(std::min(a,b) + 1);
        const uint64_t upper = uint64_t(std::max(a,b) + 1);
        return (lower << 32) | upper;
    }


}

RTetraMesh::RTetraMesh()
    : nTetras(0)
    , lastTetra(noTetra)
    , cavityStamp(0)
    , starStamp(0)
    , walkSeed(12345)
    , coplanarTolerance(0.0)
{
}

void RTetraMesh::clear()
{
    this->points.clear();
    this->tetras.clear();
    this->freeTetras.clear();
    this->vertexTetras.clear();
    this->nTetras = 0;
    this->lastTetra = noTetra;
    this->cavityMarks.clear();
    this->cavityStamp = 0;
    this->starMarks.clear();
    this->starStamp = 0;
    this->cavityTetras.clear();
    this->cavityFaces.clear();
    this->createdTetras.clear();
}

int RTetraMesh::allocateTetra()
{
    int tetraID;
    if (!this->freeTetras.empty())
    {
        tetraID = this->freeTetras.back();
        this->freeTetras.pop_back();
    }
    else
    {
        tetraID = int(this->tetras.size());
        this->tetras.push_back(Tetra());
        this->cavityMarks.push_back(0);
        this->starMarks.push_back(0);
    }
    Tetra &rTetra = this->tetras[size_t(tetraID)];
    for (int i=0;i<4;i++)
    {
        rTetra.v[i] = ghostVertex;
        rTetra.n[i] = noTetra;
    }
    rTetra.region = noRegion;
    this->nTetras++;
    return tetraID;
}

int RTetraMesh::addPoint(double x, double y, double z)
{
    Point point;
    point.x[0] = x;
    point.x[1] = y;
    point.x[2] = z;
    this->points.push_back(point);
    this->vertexTetras.push_back(noTetra);
    return int(this->points.size()) - 1;
}

int RTetraMesh::getNPoints() const
{
    return int(this->points.size());
}

int RTetraMesh::getNTetraSlots() const
{
    return int(this->tetras.size());
}

uint64_t RTetraMesh::getNTetras() const
{
    return this->nTetras;
}

void RTetraMesh::setRegion(int tetraID, int region)
{
    this->tetras[size_t(tetraID)].region = region;
}

int RTetraMesh::getVertexTetra(int pointID) const
{
    return this->vertexTetras[size_t(pointID)];
}

void RTetraMesh::createFirstTetra(int a, int b, int c, int d)
{
    const double orientation = RGeometricPredicates::orient3d(this->getPoint(a),this->getPoint(b),this->getPoint(c),this->getPoint(d));
    if (orientation == 0.0)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,"Points of the first tetrahedron are coplanar.");
    }
    if (orientation < 0.0)
    {
        std::swap(a,b);
    }

    int tetraIDs[5];

    tetraIDs[0] = this->allocateTetra();
    {
        Tetra &rTetra = this->tetras[size_t(tetraIDs[0])];
        rTetra.v[0] = a;
        rTetra.v[1] = b;
        rTetra.v[2] = c;
        rTetra.v[3] = d;
    }

    // Ghost tetrahedron on each face.
    for (int i=0;i<4;i++)
    {
        tetraIDs[i+1] = this->allocateTetra();
        const Tetra first = this->tetras[size_t(tetraIDs[0])];
        Tetra &rGhost = this->tetras[size_t(tetraIDs[i+1])];
        rGhost.v[0] = first.v[faceVertices[i][0]];
        rGhost.v[1] = first.v[faceVertices[i][2]];
        rGhost.v[2] = first.v[faceVertices[i][1]];
        rGhost.v[3] = ghostVertex;
        rGhost.n[3] = tetraIDs[0];
        this->tetras[size_t(tetraIDs[0])].n[i] = tetraIDs[i+1];
    }

    // Ghost tetrahedra are neighbors of each other.
    for (int i=1;i<5;i++)
    {
        for (int fi=0;fi<3;fi++)
        {
            const Tetra &rTetraI = this->tetras[size_t(tetraIDs[i])];
            const uint64_t keyI = findEdgeKey(rTetraI.v[(fi+1)%3],rTetraI.v[(fi+2)%3]);
            for (int j=1;j<5;j++)
            {
                if (i == j)
                {
                    continue;
                }
                const Tetra &rTetraJ = this->tetras[size_t(tetraIDs[j])];
                for (int fj=0;fj<3;fj++)
                {
                    if (keyI == findEdgeKey(rTetraJ.v[(fj+1)%3],rTetraJ.v[(fj+2)%3]))
                    {
                        this->tetras[size_t(tetraIDs[i])].n[fi] = tetraIDs[j];
                    }
                }
            }
        }
    }

    this->vertexTetras[size_t(a)] = tetraIDs[0];
    this->vertexTetras[size_t(b)] = tetraIDs[0];
    this->vertexTetras[size_t(c)] = tetraIDs[0];
    this->vertexTetras[size_t(d)] = tetraIDs[0];
    this->lastTetra = tetraIDs[0];
}

namespace
{

    //! Lexicographical order of points, which symbolic perturbation is based on.
    inline bool isPointLess(const double *p, const double *q)
    {
        if (p[0] != q[0])
        {
            return p[0] < q[0];
        }
        if (p[1] != q[1])
        {
            return p[1] < q[1];
        }
        return p[2] < q[2];
    }

    //! Return true if point p lies inside of a sphere passing through positively oriented points p0, p1, p2 and p3.
    //! If all five points lie on one sphere, they are treated as if they were moved a bit, each by a different
    //! amount given by their order (O. Devillers, M. Teillaud: Perturbations for Delaunay and weighted Delaunay
    //! 3D triangulations). The result is then the same as for points in general position, which keeps a set
    //! of tetrahedra conflicting with a point connected and without holes.
    bool isInsideSphere(const double *p0, const double *p1, const double *p2, const double *p3, const double *p)
    {
        const double side = RGeometricPredicates::inSphere(p0,p1,p2,p3,p);
        if (side != 0.0)
        {
            return side > 0.0;
        }

        const double *points[5] = {p0,p1,p2,p3,p};
        std::sort(points,points+5,isPointLess);

        for (int i=4;i>2;i--)
        {
            if (points[i] == p)
            {
                return false;
            }
            double orientation = 0.0;
            if (points[i] == p3)
            {
                orientation = RGeometricPredicates::orient3d(p0,p1,p2,p);
            }
            else if (points[i] == p2)
            {
                orientation = RGeometricPredicates::orient3d(p0,p1,p,p3);
            }
            else if (points[i] == p1)
            {
                orientation = RGeometricPredicates::orient3d(p0,p,p2,p3);
            }
            else
            {
                orientation = RGeometricPredicates::orient3d(p,p1,p2,p3);
            }
            if (orientation != 0.0)
            {
                return orientation > 0.0;
            }
        }
        return false;
    }

    //! Return true if point p lies inside of a circle passing through points p0, p1 and p2, which it is coplanar with.
    //! Point pa is any point which does not lie in their plane.
    //! Points lying on one circle are treated the same way as points lying on one sphere.
    bool isInsideCircle(const double *p0, const double *p1, const double *p2, const double *pa, const double *p)
    {
        // The circle is where a sphere passing through all four points cuts the plane.
        const double apexSide = RGeometricPredicates::orient3d(p0,p1,p2,pa);
        const double side = apexSide > 0.0 ? RGeometricPredicates::inSphere(p0,p1,p2,pa,p)
                                           : RGeometricPredicates::inSphere(p0,p2,p1,pa,p);
        if (side != 0.0)
        {
            return side > 0.0;
        }

        const double *points[4] = {p0,p1,p2,p};
        std::sort(points,points+4,isPointLess);

        // Orientation of three points in the plane is the same as their orientation with the point out of it.
        for (int i=3;i>0;i--)
        {
            if (points[i] == p)
            {
                return false;
            }
            double orientation = 0.0;
            if (points[i] == p2)
            {
                orientation = RGeometricPredicates::orient3d(p0,p1,p,pa);
            }
            else if (points[i] == p1)
            {
                orientation = RGeometricPredicates::orient3d(p0,p,p2,pa);
            }
            else
            {
                orientation = RGeometricPredicates::orient3d(p,p1,p2,pa);
            }
            if (orientation != 0.0)
            {
                return (orientation > 0.0) == (apexSide > 0.0);
            }
        }
        return false;
    }

}

void RTetraMesh::setCoplanarTolerance(double coplanarTolerance)
{
    this->coplanarTolerance = coplanarTolerance;
}

double RTetraMesh::findOrientation(const double *pa, const double *pb, const double *pc, const double *pd) const
{
    const double orientation = RGeometricPredicates::orient3d(pa,pb,pc,pd);
    if (this->coplanarTolerance > 0.0 && orientation != 0.0)
    {
        // The orientation is a product of twice the area of the triangle and of a distance of the point from it.
        const double ab[3] = {pb[0]-pa[0],pb[1]-pa[1],pb[2]-pa[2]};
        const double ac[3] = {pc[0]-pa[0],pc[1]-pa[1],pc[2]-pa[2]};
        const double n[3] = {ab[1]*ac[2]-ab[2]*ac[1],ab[2]*ac[0]-ab[0]*ac[2],ab[0]*ac[1]-ab[1]*ac[0]};
        const double n2 = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
        if (std::fabs(orientation) <= this->coplanarTolerance * std::pow(n2,0.75))
        {
            return 0.0;
        }
    }
    return orientation;
}

bool RTetraMesh::isInConflict(int tetraID, const double *p) const
{
    const Tetra &rTetra = this->tetras[size_t(tetraID)];
    const int ghostIndex = this->findGhostIndex(tetraID);
    if (ghostIndex < 0)
    {
        return isInsideSphere(this->getPoint(rTetra.v[0]),
                              this->getPoint(rTetra.v[1]),
                              this->getPoint(rTetra.v[2]),
                              this->getPoint(rTetra.v[3]),
                              p);
    }

    // Point conflicts with a ghost tetrahedron if it lies beyond its face of the convex hull.
    const int *pFace = faceVertices[ghostIndex];
    const double *p0 = this->getPoint(rTetra.v[pFace[0]]);
    const double *p1 = this->getPoint(rTetra.v[pFace[1]]);
    const double *p2 = this->getPoint(rTetra.v[pFace[2]]);
    const double orientation = this->findOrientation(p0,p1,p2,p);
    if (orientation != 0.0)
    {
        return orientation > 0.0;
    }

    // Point lies in plane of the face. Then it conflicts if it lies inside of the circumscribed circle of the face.
    const int neighborID = rTetra.n[ghostIndex];
    const Tetra &rNeighbor = this->tetras[size_t(neighborID)];
    for (int i=0;i<4;i++)
    {
        if (rNeighbor.n[i] == tetraID)
        {
            return isInsideCircle(p0,p1,p2,this->getPoint(rNeighbor.v[i]),p);
        }
    }
    return false;
}

int RTetraMesh::locate(const double *p, int startTetraID) const
{
    int tetraID = startTetraID;
    if (tetraID == noTetra || !this->isAlive(tetraID))
    {
        tetraID = this->lastTetra;
    }
    if (tetraID == noTetra)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,"Tetrahedralization is empty.");
    }
    {
        const int ghostIndex = this->findGhostIndex(tetraID);
        if (ghostIndex >= 0)
        {
            tetraID = this->tetras[size_t(tetraID)].n[ghostIndex];
        }
    }

    // Walk through faces which the point lies behind.
    int previousID = noTetra;
    const uint64_t maxSteps = this->nTetras + 16;
    for (uint64_t step=0;step<maxSteps;step++)
    {
        if (this->isGhost(tetraID))
        {
            return tetraID;
        }
        const Tetra &rTetra = this->tetras[size_t(tetraID)];
        int nextID = noTetra;
        // The first tested face changes, a walk in a mesh which is not Delaunay could go around in a circle otherwise.
        this->walkSeed = this->walkSeed * 1664525u + 1013904223u;
        const int offset = int(this->walkSeed >> 30);
        for (int k=0;k<4;k++)
        {
            const int i = (k + offset) & 3;
            if (rTetra.n[i] == previousID)
            {
                continue;
            }
            const int *pFace = faceVertices[i];
            if (RGeometricPredicates::orient3d(this->getPoint(rTetra.v[pFace[0]]),
                                               this->getPoint(rTetra.v[pFace[1]]),
                                               this->getPoint(rTetra.v[pFace[2]]),
                                               p) < 0.0)
            {
                // A point is treated as lying in plane of a face of the convex hull if it is close to it,
                // then it does not have to conflict with the ghost tetrahedron behind the face.
                if (this->isGhost(rTetra.n[i]) && !this->isInConflict(rTetra.n[i],p))
                {
                    continue;
                }
                nextID = rTetra.n[i];
                break;
            }
        }
        if (nextID == noTetra)
        {
            return tetraID;
        }
        previousID = tetraID;
        tetraID = nextID;
    }

    // The walk did not end, all tetrahedra are searched.
    int ghostID = noTetra;
    for (int i=0;i<int(this->tetras.size());i++)
    {
        if (!this->isAlive(i))
        {
            continue;
        }
        if (this->isGhost(i))
        {
            if (ghostID == noTetra && this->isInConflict(i,p))
            {
                ghostID = i;
            }
            continue;
        }
        const Tetra &rTetra = this->tetras[size_t(i)];
        bool inside = true;
        for (int j=0;j<4 && inside;j++)
        {
            const int *pFace = faceVertices[j];
            inside = RGeometricPredicates::orient3d(this->getPoint(rTetra.v[pFace[0]]),
                                                    this->getPoint(rTetra.v[pFace[1]]),
                                                    this->getPoint(rTetra.v[pFace[2]]),
                                                    p) >= 0.0;
        }
        if (inside)
        {
            return i;
        }
    }
    if (ghostID == noTetra)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,"Failed to locate a point in tetrahedralization.");
    }
    return ghostID;
}

int RTetraMesh::findDuplicateVertex(int tetraID, const double *p) const
{
    const Tetra &rTetra = this->tetras[size_t(tetraID)];
    for (int i=0;i<4;i++)
    {
        if (rTetra.v[i] == ghostVertex)
        {
            continue;
        }
        const double *q = this->getPoint(rTetra.v[i]);
        if (q[0] == p[0] && q[1] == p[1] && q[2] == p[2])
        {
            return rTetra.v[i];
        }
    }
    return -1;
}

bool RTetraMesh::findCavity(const double *p, int seedTetraID, const std::vector<int> *pForcedTetras, const std::function<bool(int,int)> *pIsBlocked)
{
    this->cavityTetras.clear();
    this->cavityFaces.clear();
    this->cavityLinks.clear();

    // Marks of tetrahedra: cavityStamp for those in the cavity, cavityStamp + 1 for those which are not
    // and cavityStamp + 2 for those which must be in the cavity.
    this->cavityStamp += 4;
    const unsigned int insideMark = this->cavityStamp;
    const unsigned int outsideMark = this->cavityStamp + 1;
    const unsigned int forcedMark = this->cavityStamp + 2;
    const unsigned int visitedMark = this->cavityStamp + 3;

    if (pForcedTetras)
    {
        for (int tetraID : *pForcedTetras)
        {
            if (this->cavityMarks[size_t(tetraID)] != forcedMark)
            {
                this->cavityMarks[size_t(tetraID)] = forcedMark;
                this->cavityTetras.push_back(tetraID);
            }
        }
    }
    // Forced tetrahedra are the first ones in the list.
    const size_t nForced = this->cavityTetras.size();

    if (this->cavityMarks[size_t(seedTetraID)] != forcedMark)
    {
        if (this->isInConflict(seedTetraID,p))
        {
            this->cavityMarks[size_t(seedTetraID)] = insideMark;
            this->cavityTetras.push_back(seedTetraID);
        }
        else if (this->cavityTetras.empty())
        {
            return false;
        }
    }

    // Tetrahedra conflicting with the point, which can be reached without crossing a blocked face.
    for (size_t i=0;i<this->cavityTetras.size();i++)
    {
        const int tetraID = this->cavityTetras[i];
        for (int j=0;j<4;j++)
        {
            const int neighborID = this->tetras[size_t(tetraID)].n[j];
            if (this->isInCavity(neighborID) || this->cavityMarks[size_t(neighborID)] == outsideMark)
            {
                continue;
            }
            if (pIsBlocked && (*pIsBlocked)(tetraID,j))
            {
                continue;
            }
            if (this->isInConflict(neighborID,p))
            {
                this->cavityMarks[size_t(neighborID)] = insideMark;
                this->cavityTetras.push_back(neighborID);
            }
            else
            {
                this->cavityMarks[size_t(neighborID)] = outsideMark;
            }
        }
    }

    // Tetrahedra are removed from the cavity until the point can see all faces of its boundary
    // and until no blocked face is inside of it.
    std::vector<int> stack;
    std::vector<int> component;
    std::vector<int> bestComponent;

    while (true)
    {
        bool changed = false;

        for (int tetraID : this->cavityTetras)
        {
            if (!this->isInCavity(tetraID))
            {
                continue;
            }
            const Tetra &rTetra = this->tetras[size_t(tetraID)];
            for (int j=0;j<4;j++)
            {
                const int neighborID = rTetra.n[j];
                const int *pFace = faceVertices[j];
                const int a = rTetra.v[pFace[0]];
                const int b = rTetra.v[pFace[1]];
                const int c = rTetra.v[pFace[2]];
                const bool finite = a != ghostVertex && b != ghostVertex && c != ghostVertex;

                int removeID = noTetra;
                if (this->isInCavity(neighborID))
                {
                    if (!pIsBlocked || !(*pIsBlocked)(tetraID,j))
                    {
                        continue;
                    }
                    // Tetrahedron on the side of the face where the point is stays.
                    const double side = finite ? RGeometricPredicates::orient3d(this->getPoint(a),this->getPoint(b),this->getPoint(c),p) : 1.0;
                    removeID = side < 0.0 ? tetraID : neighborID;
                    if (this->cavityMarks[size_t(removeID)] == forcedMark)
                    {
                        removeID = removeID == tetraID ? neighborID : tetraID;
                    }
                }
                else if (finite && !(this->findOrientation(this->getPoint(a),this->getPoint(b),this->getPoint(c),p) > 0.0))
                {
                    removeID = tetraID;
                }
                if (removeID == noTetra)
                {
                    continue;
                }
                if (this->cavityMarks[size_t(removeID)] == forcedMark)
                {
                    this->cavityTetras.clear();
                    return false;
                }
                this->cavityMarks[size_t(removeID)] = outsideMark;
                changed = true;
                if (removeID == tetraID)
                {
                    break;
                }
            }
        }

        if (!changed)
        {
            break;
        }

        // Only one group of connected tetrahedra is kept: the one with forced tetrahedra,
        // otherwise the one with the seed and if the seed was removed, then the largest one.
        bestComponent.clear();
        int bestRank = -1;
        int nForcedGroups = 0;
        for (int firstID : this->cavityTetras)
        {
            if (!this->isInCavity(firstID))
            {
                continue;
            }
            component.clear();
            stack.clear();
            bool hasForced = false;
            bool hasSeed = false;

            auto visit = [&](int tetraID)
            {
                hasForced = hasForced || this->cavityMarks[size_t(tetraID)] == forcedMark;
                hasSeed = hasSeed || tetraID == seedTetraID;
                this->cavityMarks[size_t(tetraID)] = visitedMark;
                component.push_back(tetraID);
                stack.push_back(tetraID);
            };

            visit(firstID);
            while (!stack.empty())
            {
                const int tetraID = stack.back();
                stack.pop_back();
                for (int j=0;j<4;j++)
                {
                    const int neighborID = this->tetras[size_t(tetraID)].n[j];
                    if (this->isInCavity(neighborID))
                    {
                        visit(neighborID);
                    }
                }
            }
            if (hasForced)
            {
                nForcedGroups++;
            }
            const int rank = hasForced ? 2 : (hasSeed ? 1 : 0);
            if (rank > bestRank || (rank == bestRank && component.size() > bestComponent.size()))
            {
                bestRank = rank;
                bestComponent = component;
            }
        }
        if (nForcedGroups > 1 || bestComponent.empty())
        {
            this->cavityTetras.clear();
            return false;
        }
        for (int tetraID : this->cavityTetras)
        {
            if (this->cavityMarks[size_t(tetraID)] == visitedMark)
            {
                this->cavityMarks[size_t(tetraID)] = outsideMark;
            }
        }
        for (int tetraID : bestComponent)
        {
            this->cavityMarks[size_t(tetraID)] = insideMark;
        }
        // All forced tetrahedra are in the kept group.
        for (size_t i=0;i<nForced;i++)
        {
            this->cavityMarks[size_t(this->cavityTetras[i])] = forcedMark;
        }
    }

    size_t nKept = 0;
    for (size_t i=0;i<this->cavityTetras.size();i++)
    {
        if (this->isInCavity(this->cavityTetras[i]))
        {
            this->cavityTetras[nKept++] = this->cavityTetras[i];
        }
    }
    this->cavityTetras.resize(nKept);
    if (this->cavityTetras.empty())
    {
        return false;
    }

    for (int tetraID : this->cavityTetras)
    {
        for (int j=0;j<4;j++)
        {
            if (!this->isInCavity(this->tetras[size_t(tetraID)].n[j]))
            {
                this->cavityFaces.push_back(Face{tetraID,j});
            }
        }
    }

    // Each edge of the boundary must be shared by two of its faces, so that new tetrahedra can be connected.
    this->cavityLinks.reserve(3 * this->cavityFaces.size());
    for (size_t i=0;i<this->cavityFaces.size();i++)
    {
        const Tetra &rTetra = this->tetras[size_t(this->cavityFaces[i].tetra)];
        const int *pFace = faceVertices[this->cavityFaces[i].face];
        const int a = rTetra.v[pFace[0]];
        const int b = rTetra.v[pFace[1]];
        const int c = rTetra.v[pFace[2]];
        this->cavityLinks.push_back(FaceLink{findEdgeKey(b,c),int(i),0});
        this->cavityLinks.push_back(FaceLink{findEdgeKey(a,c),int(i),1});
        this->cavityLinks.push_back(FaceLink{findEdgeKey(a,b),int(i),2});
    }
    std::sort(this->cavityLinks.begin(),this->cavityLinks.end(),[](const FaceLink &l1, const FaceLink &l2) { return l1.key < l2.key; });
    for (size_t i=0;i<this->cavityLinks.size();i+=2)
    {
        if (i + 1 >= this->cavityLinks.size()
            || this->cavityLinks[i].key != this->cavityLinks[i+1].key
            || (i + 2 < this->cavityLinks.size() && this->cavityLinks[i+2].key == this->cavityLinks[i].key))
        {
            this->cavityTetras.clear();
            this->cavityFaces.clear();
            this->cavityLinks.clear();
            return false;
        }
    }
    return true;
}

const std::vector<int> &RTetraMesh::getCavityTetras() const
{
    return this->cavityTetras;
}

const std::vector<RTetraMesh::Face> &RTetraMesh::getCavityFaces() const
{
    return this->cavityFaces;
}

void RTetraMesh::fillCavity(int pointID)
{
    this->createdTetras.clear();

    for (const Face &rFace : this->cavityFaces)
    {
        const Tetra oldTetra = this->tetras[size_t(rFace.tetra)];
        const int *pFace = faceVertices[rFace.face];
        const int outsideID = oldTetra.n[rFace.face];

        const int newID = this->allocateTetra();
        Tetra &rNew = this->tetras[size_t(newID)];
        rNew.v[0] = oldTetra.v[pFace[0]];
        rNew.v[1] = oldTetra.v[pFace[1]];
        rNew.v[2] = oldTetra.v[pFace[2]];
        rNew.v[3] = pointID;
        rNew.n[3] = outsideID;
        rNew.region = oldTetra.region;

        Tetra &rOutside = this->tetras[size_t(outsideID)];
        for (int i=0;i<4;i++)
        {
            if (rOutside.n[i] == rFace.tetra)
            {
                rOutside.n[i] = newID;
                break;
            }
        }

        for (int i=0;i<3;i++)
        {
            if (rNew.v[i] != ghostVertex)
            {
                this->vertexTetras[size_t(rNew.v[i])] = newID;
            }
        }
        this->vertexTetras[size_t(pointID)] = newID;
        this->createdTetras.push_back(newID);
    }

    // New tetrahedra sharing an edge of the cavity boundary are neighbors.
    for (size_t i=0;i+1<this->cavityLinks.size();i+=2)
    {
        const FaceLink &rLink1 = this->cavityLinks[i];
        const FaceLink &rLink2 = this->cavityLinks[i+1];
        const int tetraID1 = this->createdTetras[size_t(rLink1.face)];
        const int tetraID2 = this->createdTetras[size_t(rLink2.face)];
        this->tetras[size_t(tetraID1)].n[rLink1.side] = tetraID2;
        this->tetras[size_t(tetraID2)].n[rLink2.side] = tetraID1;
    }

    for (int tetraID : this->cavityTetras)
    {
        this->tetras[size_t(tetraID)].v[0] = deadVertex;
        this->freeTetras.push_back(tetraID);
        this->nTetras--;
    }
    this->cavityTetras.clear();
    this->cavityFaces.clear();
    this->cavityLinks.clear();
    // Marks of the cavity must not be found on tetrahedra which take positions of the removed ones.
    this->cavityStamp += 4;

    for (int tetraID : this->createdTetras)
    {
        if (!this->isGhost(tetraID))
        {
            this->lastTetra = tetraID;
            this->vertexTetras[size_t(pointID)] = tetraID;
            break;
        }
    }
}

const std::vector<int> &RTetraMesh::getCreatedTetras() const
{
    return this->createdTetras;
}

bool RTetraMesh::insertVertex(int pointID, int hintTetraID, int &duplicateID)
{
    const double *p = this->getPoint(pointID);
    const int seedID = this->locate(p,hintTetraID);

    duplicateID = this->findDuplicateVertex(seedID,p);
    if (duplicateID >= 0)
    {
        return false;
    }
    if (!this->findCavity(p,seedID))
    {
        // The point is so close to a vertex, that it can not see faces around it at all.
        // It is treated as a duplicate of the closest vertex of the tetrahedron it lies in.
        const Tetra &rTetra = this->tetras[size_t(seedID)];
        double minDistance = 0.0;
        for (int i=0;i<4;i++)
        {
            if (rTetra.v[i] == ghostVertex)
            {
                continue;
            }
            const double *q = this->getPoint(rTetra.v[i]);
            const double distance = std::hypot(q[0]-p[0],q[1]-p[1],q[2]-p[2]);
            if (duplicateID < 0 || distance < minDistance)
            {
                minDistance = distance;
                duplicateID = rTetra.v[i];
            }
        }
        return false;
    }
    this->fillCavity(pointID);
    return true;
}

void RTetraMesh::findStar(int pointID, std::vector<int> &star) const
{
    star.clear();

    const int firstID = this->vertexTetras[size_t(pointID)];
    if (firstID == noTetra)
    {
        return;
    }

    const unsigned int mark = ++this->starStamp;
    this->starMarks[size_t(firstID)] = mark;
    star.push_back(firstID);

    for (size_t i=0;i<star.size();i++)
    {
        const Tetra &rTetra = this->tetras[size_t(star[i])];
        for (int j=0;j<4;j++)
        {
            // Neighbor opposite to the vertex does not have it.
            if (rTetra.v[j] == pointID)
            {
                continue;
            }
            const int neighborID = rTetra.n[j];
            if (this->starMarks[size_t(neighborID)] != mark)
            {
                this->starMarks[size_t(neighborID)] = mark;
                star.push_back(neighborID);
            }
        }
    }
}

int RTetraMesh::findEdge(int a, int b) const
{
    this->findStar(a,this->starBuffer);
    for (int tetraID : this->starBuffer)
    {
        const Tetra &rTetra = this->tetras[size_t(tetraID)];
        if (rTetra.v[0] == b || rTetra.v[1] == b || rTetra.v[2] == b || rTetra.v[3] == b)
        {
            return tetraID;
        }
    }
    return noTetra;
}

bool RTetraMesh::findFace(int a, int b, int c, Face &face) const
{
    this->findStar(a,this->starBuffer);
    for (int tetraID : this->starBuffer)
    {
        const Tetra &rTetra = this->tetras[size_t(tetraID)];
        int opposite = -1;
        int nFound = 0;
        for (int i=0;i<4;i++)
        {
            if (rTetra.v[i] == a || rTetra.v[i] == b || rTetra.v[i] == c)
            {
                nFound++;
            }
            else
            {
                opposite = i;
            }
        }
        if (nFound == 3)
        {
            face.tetra = tetraID;
            face.face = opposite;
            return true;
        }
    }
    return false;
}
