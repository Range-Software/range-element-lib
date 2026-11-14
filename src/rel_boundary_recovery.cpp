#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>

#include <rbl_error.h>

#include "rel_boundary_recovery.h"
#include "rel_tetra_mesh_quality.h"

namespace
{

    //! Relative distance below which a point is considered to lie on an edge of a subface.
    const double onEdgeTolerance = 1.0e-9;
    //! Relative margin by which a vertex has to lie inside of a circle to flip an edge.
    const double flipTolerance = 1.0e-9;
    //! Relative distance from a circle below which a vertex is considered to lie on it.
    const double cocircularTolerance = 1.0e-6;

    inline uint64_t findEdgeKey(int a, int b)
    {
        return (uint64_t(std::min(a,b)) << 32) | uint64_t(std::max(a,b));
    }

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

    inline double length(const double *a)
    {
        return std::sqrt(dot(a,a));
    }

    inline double distance(const double *a, const double *b)
    {
        double r[3];
        subtract(a,b,r);
        return length(r);
    }

    //! Normal of a triangle, its length is twice the area of the triangle.
    inline void findNormal(const double *a, const double *b, const double *c, double *n)
    {
        double ab[3], ac[3];
        subtract(b,a,ab);
        subtract(c,a,ac);
        cross(ab,ac,n);
    }

    //! Center and radius of a circumscribed circle of a triangle. Return false if the triangle is degenerated.
    bool findCircumcircle(const double *a, const double *b, const double *c, double *center, double &radius)
    {
        double ab[3], ac[3], n[3], nab[3], acn[3];
        subtract(b,a,ab);
        subtract(c,a,ac);
        cross(ab,ac,n);
        const double n2 = dot(n,n);
        if (n2 <= 0.0)
        {
            return false;
        }
        cross(n,ab,nab);
        cross(ac,n,acn);
        const double ab2 = dot(ab,ab);
        const double ac2 = dot(ac,ac);
        double r[3];
        for (int i=0;i<3;i++)
        {
            r[i] = (ac2 * nab[i] + ab2 * acn[i]) / (2.0 * n2);
            center[i] = a[i] + r[i];
        }
        radius = length(r);
        return std::isfinite(radius);
    }

    //! Return true if point lies inside of a sphere which the segment is a diameter of.
    inline bool isInDiametralSphere(const double *p, const double *a, const double *b)
    {
        double pa[3], pb[3];
        subtract(a,p,pa);
        subtract(b,p,pb);
        return dot(pa,pb) < 0.0;
    }

}

RBoundaryRecovery::RBoundaryRecovery(RTetraMesh &mesh, double tolerance)
    : mesh(mesh)
    , tolerance(std::max(tolerance,1.0e-12))
    , minLength(0.0)
    , maxSteinerVertices(std::numeric_limits<uint64_t>::max())
    , vertexTypes(size_t(mesh.getNPoints()),VertexInput)
    , boundaryVertices(size_t(mesh.getNPoints()),0)
    , vertexSizes(size_t(mesh.getNPoints()),0.0)
    , insertionRadii(size_t(mesh.getNPoints()),0.0)
    , vertexParents(size_t(mesh.getNPoints()))
    , replacedStamp(0)
    , facetStamp(0)
    , collectCreatedTetras(false)
{
    for (int i=0;i<mesh.getNPoints();i++) this->vertexParents[size_t(i)]={i};
    this->blockedFunction = [this](int tetraID, int face) -> bool
    {
        return this->isFaceBlocked(tetraID,face);
    };
    for (uint64_t &rCount : this->nSteinerVertices)
    {
        rCount = 0;
    }
}

void RBoundaryRecovery::setMinLength(double minLength)
{
    this->minLength = minLength;
}

void RBoundaryRecovery::setMaxSteinerVertices(uint64_t maxSteinerVertices)
{
    this->maxSteinerVertices = maxSteinerVertices;
}

uint64_t RBoundaryRecovery::getNSteinerVertices(VertexType type) const
{
    return this->nSteinerVertices[type];
}

uint64_t RBoundaryRecovery::getNSteinerVertices() const
{
    uint64_t count = 0;
    for (uint64_t typeCount : this->nSteinerVertices)
    {
        count += typeCount;
    }
    return count;
}

double RBoundaryRecovery::getVertexSize(int vertexID) const
{
    return this->vertexSizes[size_t(vertexID)];
}

void RBoundaryRecovery::setVertexSize(int vertexID, double size)
{
    this->vertexSizes[size_t(vertexID)] = size;
}

RBoundaryRecovery::VertexType RBoundaryRecovery::getVertexType(int vertexID) const
{
    return VertexType(this->vertexTypes[size_t(vertexID)]);
}

int RBoundaryRecovery::addSubSegment(int a, int b, int marker)
{
    if (a == b)
    {
        return -1;
    }
    const uint64_t key = findEdgeKey(a,b);
    auto iter = this->segmentMap.find(key);
    if (iter != this->segmentMap.end())
    {
        if (marker != noMarker)
        {
            this->segments[size_t(iter->second)].marker = marker;
        }
        return iter->second;
    }

    SubSegment segment;
    segment.v[0] = a;
    segment.v[1] = b;
    segment.marker = marker;
    segment.alive = true;
    segment.queued = false;
    this->segments.push_back(segment);

    const int segmentID = int(this->segments.size()) - 1;
    this->segmentMap[key] = segmentID;
    this->boundaryVertices[size_t(a)] = 1;
    this->boundaryVertices[size_t(b)] = 1;
    return segmentID;
}

void RBoundaryRecovery::removeSubSegment(int segmentID)
{
    SubSegment &rSegment = this->segments[size_t(segmentID)];
    rSegment.alive = false;
    this->segmentMap.erase(findEdgeKey(rSegment.v[0],rSegment.v[1]));
}

int RBoundaryRecovery::addSubFace(int a, int b, int c, int marker)
{
    if (a == b || b == c || a == c)
    {
        return -1;
    }
    FaceKey key = {{a,b,c}};
    std::sort(key.v,key.v+3);
    if (this->faceMap.find(key) != this->faceMap.end())
    {
        return -1;
    }

    SubFace face;
    face.v[0] = a;
    face.v[1] = b;
    face.v[2] = c;
    face.marker = marker;
    face.alive = true;
    face.queued = false;
    this->faces.push_back(face);
    this->replacedMarks.push_back(0);
    this->facetMarks.push_back(0);
    this->facetMatchPoints.push_back(-1);

    const int faceID = int(this->faces.size()) - 1;
    this->faceMap[key] = faceID;
    for (int i=0;i<3;i++)
    {
        this->edgeFaces[findEdgeKey(face.v[i],face.v[(i+1)%3])].push_back(faceID);
        this->boundaryVertices[size_t(face.v[i])] = 1;
    }
    return faceID;
}

void RBoundaryRecovery::removeSubFace(int faceID)
{
    SubFace &rFace = this->faces[size_t(faceID)];
    rFace.alive = false;

    FaceKey key = {{rFace.v[0],rFace.v[1],rFace.v[2]}};
    std::sort(key.v,key.v+3);
    this->faceMap.erase(key);

    for (int i=0;i<3;i++)
    {
        auto iter = this->edgeFaces.find(findEdgeKey(rFace.v[i],rFace.v[(i+1)%3]));
        if (iter == this->edgeFaces.end())
        {
            continue;
        }
        std::vector<int> &rFaceIDs = iter->second;
        rFaceIDs.erase(std::remove(rFaceIDs.begin(),rFaceIDs.end(),faceID),rFaceIDs.end());
        if (rFaceIDs.empty())
        {
            this->edgeFaces.erase(iter);
        }
    }
}

void RBoundaryRecovery::queueSubSegment(int segmentID)
{
    if (segmentID < 0)
    {
        return;
    }
    SubSegment &rSegment = this->segments[size_t(segmentID)];
    if (rSegment.alive && !rSegment.queued)
    {
        rSegment.queued = true;
        this->segmentQueue.push_back(segmentID);
    }
}

void RBoundaryRecovery::queueSubFace(int faceID)
{
    if (faceID < 0)
    {
        return;
    }
    SubFace &rFace = this->faces[size_t(faceID)];
    if (rFace.alive && !rFace.queued)
    {
        rFace.queued = true;
        this->faceQueue.push_back(faceID);
    }
}

void RBoundaryRecovery::queueCavity()
{
    for (int tetraID : this->mesh.getCavityTetras())
    {
        const RTetraMesh::Tetra &rTetra = this->mesh.getTetra(tetraID);
        bool boundary[4];
        int nBoundary = 0;
        for (int i=0;i<4;i++)
        {
            boundary[i] = rTetra.v[i] != RTetraMesh::ghostVertex && this->boundaryVertices[size_t(rTetra.v[i])];
            if (boundary[i])
            {
                nBoundary++;
            }
        }
        if (nBoundary < 2)
        {
            continue;
        }
        for (int i=0;i<6;i++)
        {
            const int i1 = RTetraMesh::edgeVertices[i][0];
            const int i2 = RTetraMesh::edgeVertices[i][1];
            if (boundary[i1] && boundary[i2])
            {
                this->queueSubSegment(this->findSubSegment(rTetra.v[i1],rTetra.v[i2]));
            }
        }
        if (nBoundary < 3)
        {
            continue;
        }
        for (int i=0;i<4;i++)
        {
            const int *pFace = RTetraMesh::faceVertices[i];
            if (boundary[pFace[0]] && boundary[pFace[1]] && boundary[pFace[2]])
            {
                this->queueSubFace(this->findSubFace(rTetra.v[pFace[0]],rTetra.v[pFace[1]],rTetra.v[pFace[2]]));
            }
        }
    }
}

double RBoundaryRecovery::interpolateSize(const double *p, int tetraID) const
{
    const RTetraMesh::Tetra &rTetra = this->mesh.getTetra(tetraID);
    double weightSum = 0.0;
    double sizeSum = 0.0;
    for (int i=0;i<4;i++)
    {
        if (rTetra.v[i] == RTetraMesh::ghostVertex)
        {
            continue;
        }
        const double size = this->vertexSizes[size_t(rTetra.v[i])];
        if (size <= 0.0)
        {
            continue;
        }
        const double d = distance(p,this->mesh.getPoint(rTetra.v[i]));
        if (d <= 0.0)
        {
            return size;
        }
        weightSum += 1.0 / d;
        sizeSum += size / d;
    }
    return weightSum > 0.0 ? sizeSum / weightSum : 0.0;
}

bool RBoundaryRecovery::isFaceBlocked(int tetraID, int face) const
{
    const int faceID = this->findSubFace(RTetraMesh::Face{tetraID,face});
    return faceID >= 0 && this->replacedMarks[size_t(faceID)] != this->replacedStamp;
}

void RBoundaryRecovery::findFacetCavity(const double *p, const std::vector<int> &startFaceIDs, uint64_t exceptKey, std::vector<int> &removedFaceIDs)
{
    removedFaceIDs = startFaceIDs;

    const unsigned int mark = ++this->facetStamp;
    for (int faceID : removedFaceIDs)
    {
        this->facetMarks[size_t(faceID)] = mark;
    }

    // Neighbor on the other side of an edge, -1 if the edge is a subsegment or if there is none.
    auto findNeighbor = [this](int faceID, int a, int b) -> int
    {
        if (this->findSubSegment(a,b) >= 0)
        {
            return -1;
        }
        auto iter = this->edgeFaces.find(findEdgeKey(a,b));
        if (iter == this->edgeFaces.end() || iter->second.size() != 2)
        {
            return -1;
        }
        return iter->second[0] == faceID ? iter->second[1] : iter->second[0];
    };

    // Subfaces of the facet which circumscribed circles contain the point.
    for (size_t i=0;i<removedFaceIDs.size();i++)
    {
        const SubFace face = this->faces[size_t(removedFaceIDs[i])];
        for (int j=0;j<3;j++)
        {
            const int a = face.v[j];
            const int b = face.v[(j+1)%3];
            if (findEdgeKey(a,b) == exceptKey)
            {
                continue;
            }
            const int neighborID = findNeighbor(removedFaceIDs[i],a,b);
            if (neighborID < 0 || this->facetMarks[size_t(neighborID)] == mark)
            {
                continue;
            }
            const SubFace &rNeighbor = this->faces[size_t(neighborID)];
            double center[3], radius;
            if (!findCircumcircle(this->mesh.getPoint(rNeighbor.v[0]),this->mesh.getPoint(rNeighbor.v[1]),this->mesh.getPoint(rNeighbor.v[2]),center,radius))
            {
                continue;
            }
            if (distance(p,center) < radius * (1.0 - flipTolerance))
            {
                this->facetMarks[size_t(neighborID)] = mark;
                removedFaceIDs.push_back(neighborID);
            }
        }
    }

    // Subfaces are taken out until the point can see all edges around them, and until they are all connected.
    const size_t nStart = startFaceIDs.size();
    std::vector<int> stack;
    while (true)
    {
        bool changed = false;
        for (size_t i=0;i<removedFaceIDs.size();i++)
        {
            const int faceID = removedFaceIDs[i];
            if (this->facetMarks[size_t(faceID)] != mark)
            {
                continue;
            }
            const SubFace &rFace = this->faces[size_t(faceID)];
            const double *pv[3] = {this->mesh.getPoint(rFace.v[0]),this->mesh.getPoint(rFace.v[1]),this->mesh.getPoint(rFace.v[2])};
            double n[3];
            findNormal(pv[0],pv[1],pv[2],n);
            const double nLength = length(n);
            const double scale = std::max(distance(pv[0],pv[1]),std::max(distance(pv[1],pv[2]),distance(pv[2],pv[0])));
            for (int j=0;j<3;j++)
            {
                const int a = rFace.v[j];
                const int b = rFace.v[(j+1)%3];
                if (findEdgeKey(a,b) == exceptKey)
                {
                    continue;
                }
                const int neighborID = findNeighbor(faceID,a,b);
                if (neighborID >= 0 && this->facetMarks[size_t(neighborID)] == mark)
                {
                    continue;
                }
                double edge[3], toPoint[3], product[3];
                subtract(pv[(j+1)%3],pv[j],edge);
                subtract(p,pv[j],toPoint);
                cross(edge,toPoint,product);
                const double edgeLength = length(edge);
                if (!(nLength > 0.0) || !(dot(product,n) > onEdgeTolerance * scale * nLength * edgeLength))
                {
                    if (i < nStart)
                    {
                        // A subface which the point lies in can not be taken out.
                        removedFaceIDs.clear();
                        return;
                    }
                    this->facetMarks[size_t(faceID)] = 0;
                    changed = true;
                    break;
                }
            }
        }
        if (!changed)
        {
            break;
        }

        const unsigned int connectedMark = ++this->facetStamp;
        stack.clear();
        for (size_t i=0;i<nStart;i++)
        {
            this->facetMarks[size_t(removedFaceIDs[i])] = connectedMark;
            stack.push_back(removedFaceIDs[i]);
        }
        while (!stack.empty())
        {
            const int faceID = stack.back();
            stack.pop_back();
            const SubFace &rFace = this->faces[size_t(faceID)];
            for (int j=0;j<3;j++)
            {
                const int a = rFace.v[j];
                const int b = rFace.v[(j+1)%3];
                if (findEdgeKey(a,b) == exceptKey)
                {
                    continue;
                }
                const int neighborID = findNeighbor(faceID,a,b);
                if (neighborID >= 0 && this->facetMarks[size_t(neighborID)] == mark)
                {
                    this->facetMarks[size_t(neighborID)] = connectedMark;
                    stack.push_back(neighborID);
                }
            }
        }
        // Marks of connected subfaces are set back.
        for (int faceID : removedFaceIDs)
        {
            unsigned int &rMark = this->facetMarks[size_t(faceID)];
            rMark = rMark == connectedMark ? mark : 0;
        }
    }

    size_t nKept = 0;
    for (size_t i=0;i<removedFaceIDs.size();i++)
    {
        if (this->facetMarks[size_t(removedFaceIDs[i])] == mark)
        {
            removedFaceIDs[nKept++] = removedFaceIDs[i];
        }
    }
    removedFaceIDs.resize(nKept);
}

std::vector<RBoundaryRecovery::LoopEdge> RBoundaryRecovery::findLoopEdges(const std::vector<int> &removedFaceIDs, uint64_t exceptKey)
{
    std::vector<LoopEdge> loopEdges;

    const unsigned int mark = ++this->facetStamp;
    for (int faceID : removedFaceIDs)
    {
        this->facetMarks[size_t(faceID)] = mark;
    }
    for (int faceID : removedFaceIDs)
    {
        const SubFace &rFace = this->faces[size_t(faceID)];
        for (int j=0;j<3;j++)
        {
            const int a = rFace.v[j];
            const int b = rFace.v[(j+1)%3];
            const uint64_t key = findEdgeKey(a,b);
            if (key == exceptKey)
            {
                continue;
            }
            bool inside = false;
            if (this->findSubSegment(a,b) < 0)
            {
                auto iter = this->edgeFaces.find(key);
                if (iter != this->edgeFaces.end() && iter->second.size() == 2)
                {
                    const int neighborID = iter->second[0] == faceID ? iter->second[1] : iter->second[0];
                    inside = this->facetMarks[size_t(neighborID)] == mark;
                }
            }
            if (!inside)
            {
                loopEdges.push_back(LoopEdge{a,b,rFace.marker});
            }
        }
    }
    return loopEdges;
}

bool RBoundaryRecovery::findBoundaryCavity(const double *p, const std::vector<int> &removedFaceIDs, int a, int b, int locatedTetraID, bool exactVisibility)
{
    // Removed subfaces do not block the cavity and tetrahedra on both of their sides are a part of it.
    this->replacedStamp++;
    std::vector<int> forcedTetras;
    for (int faceID : removedFaceIDs)
    {
        this->replacedMarks[size_t(faceID)] = this->replacedStamp;
        const SubFace &rFace = this->faces[size_t(faceID)];
        RTetraMesh::Face meshFace;
        if (this->mesh.findFace(rFace.v[0],rFace.v[1],rFace.v[2],meshFace))
        {
            forcedTetras.push_back(meshFace.tetra);
            forcedTetras.push_back(this->mesh.getTetra(meshFace.tetra).n[meshFace.face]);
        }
    }
    // So are all tetrahedra around a split edge.
    if (a >= 0 && b >= 0)
    {
        std::vector<int> star;
        this->mesh.findStar(a,star);
        for (int tetraID : star)
        {
            const RTetraMesh::Tetra &rTetra = this->mesh.getTetra(tetraID);
            if (rTetra.v[0] == b || rTetra.v[1] == b || rTetra.v[2] == b || rTetra.v[3] == b)
            {
                forcedTetras.push_back(tetraID);
            }
        }
    }
    const int seedID = forcedTetras.empty() ? locatedTetraID : forcedTetras.front();
    return this->mesh.findCavity(p,seedID,&forcedTetras,&this->blockedFunction,exactVisibility);
}

bool RBoundaryRecovery::prepareVolumeCavity(const double *p, int seedTetraID)
{
    // No subface is replaced, all of them block the cavity.
    this->replacedStamp++;
    return this->mesh.findCavity(p,seedTetraID,nullptr,&this->blockedFunction);
}

int RBoundaryRecovery::findLostSubSegment(int exceptSegmentID) const
{
    // Edges which remain are those on boundary of the cavity.
    std::vector<uint64_t> boundaryKeys;
    for (const RTetraMesh::Face &rFace : this->mesh.getCavityFaces())
    {
        const RTetraMesh::Tetra &rTetra = this->mesh.getTetra(rFace.tetra);
        const int *pFace = RTetraMesh::faceVertices[rFace.face];
        for (int i=0;i<3;i++)
        {
            const int a = rTetra.v[pFace[i]];
            const int b = rTetra.v[pFace[(i+1)%3]];
            if (a != RTetraMesh::ghostVertex && b != RTetraMesh::ghostVertex
                && this->boundaryVertices[size_t(a)] && this->boundaryVertices[size_t(b)])
            {
                boundaryKeys.push_back(findEdgeKey(a,b));
            }
        }
    }
    std::sort(boundaryKeys.begin(),boundaryKeys.end());

    for (int tetraID : this->mesh.getCavityTetras())
    {
        const RTetraMesh::Tetra &rTetra = this->mesh.getTetra(tetraID);
        for (int i=0;i<6;i++)
        {
            const int a = rTetra.v[RTetraMesh::edgeVertices[i][0]];
            const int b = rTetra.v[RTetraMesh::edgeVertices[i][1]];
            if (a == RTetraMesh::ghostVertex || b == RTetraMesh::ghostVertex
                || !this->boundaryVertices[size_t(a)] || !this->boundaryVertices[size_t(b)])
            {
                continue;
            }
            const int segmentID = this->findSubSegment(a,b);
            if (segmentID < 0 || segmentID == exceptSegmentID)
            {
                continue;
            }
            if (!std::binary_search(boundaryKeys.begin(),boundaryKeys.end(),findEdgeKey(a,b)))
            {
                return segmentID;
            }
        }
    }
    return -1;
}

double RBoundaryRecovery::findCavityRadius(const double *p) const
{
    double radius = std::numeric_limits<double>::max();
    for (int tetraID : this->mesh.getCavityTetras())
    {
        const RTetraMesh::Tetra &rTetra = this->mesh.getTetra(tetraID);
        for (int i=0;i<4;i++)
        {
            if (rTetra.v[i] != RTetraMesh::ghostVertex)
            {
                radius = std::min(radius,distance(p,this->mesh.getPoint(rTetra.v[i])));
            }
        }
    }
    return radius;
}

int RBoundaryRecovery::insertCavityPoint(const double *p, VertexType type, double radius, const std::vector<int> &parents)
{
    if (this->getNSteinerVertices() >= this->maxSteinerVertices)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Maximum number of added nodes (%llu) was reached.",
                     static_cast<unsigned long long>(this->maxSteinerVertices));
    }

    const double size = this->interpolateSize(p,this->mesh.getCavityTetras().front());

    const double spacing=this->findCavityRadius(p);
    std::vector<int> ancestry;
    for (int parent : parents)
        ancestry.insert(ancestry.end(),this->vertexParents[size_t(parent)].begin(),this->vertexParents[size_t(parent)].end());
    std::sort(ancestry.begin(),ancestry.end());
    ancestry.erase(std::unique(ancestry.begin(),ancestry.end()),ancestry.end());
    this->queueCavity();

    const int vertexID = this->mesh.addPoint(p[0],p[1],p[2]);
    this->vertexTypes.push_back(static_cast<unsigned char>(type));
    this->boundaryVertices.push_back(0);
    this->vertexSizes.push_back(size);
    this->insertionRadii.push_back(radius>0.0?std::min(radius,spacing):0.0);
    this->vertexParents.push_back(ancestry);

    this->mesh.fillCavity(vertexID);

    if (this->collectCreatedTetras)
    {
        const std::vector<int> &rCreated = this->mesh.getCreatedTetras();
        this->createdTetras.insert(this->createdTetras.end(),rCreated.begin(),rCreated.end());
    }
    this->nSteinerVertices[type]++;

    return vertexID;
}

void RBoundaryRecovery::addTriangle(int a, int b, int c, int marker)
{
    if (a == b || b == c || a == c)
    {
        return;
    }
    double n[3];
    findNormal(this->mesh.getPoint(a),this->mesh.getPoint(b),this->mesh.getPoint(c),n);
    if (!(dot(n,n) > 0.0))
    {
        // Triangle without an area can not be a part of a boundary.
        return;
    }
    this->addSubFace(a,b,c,marker);
}

void RBoundaryRecovery::addEdge(int a, int b, int marker)
{
    this->addSubSegment(a,b,marker);
}

void RBoundaryRecovery::initialize()
{
    // Edges which facets are bounded by are subsegments.
    std::vector<uint64_t> facetEdges;
    for (const auto &rItem : this->edgeFaces)
    {
        const uint64_t key = rItem.first;
        const std::vector<int> &rFaceIDs = rItem.second;
        const int a = int(key >> 32);
        const int b = int(key & 0xffffffffULL);

        bool isSegment = true;
        if (rFaceIDs.size() == 2)
        {
            const SubFace &rFace1 = this->faces[size_t(rFaceIDs[0])];
            const SubFace &rFace2 = this->faces[size_t(rFaceIDs[1])];
            if (rFace1.marker == rFace2.marker)
            {
                int c = -1, d = -1;
                for (int i=0;i<3;i++)
                {
                    if (rFace1.v[i] != a && rFace1.v[i] != b)
                    {
                        c = rFace1.v[i];
                    }
                    if (rFace2.v[i] != a && rFace2.v[i] != b)
                    {
                        d = rFace2.v[i];
                    }
                }
                // Subfaces lie in one plane if their normals, taken as if both subfaces
                // went around the edge the same way, point against each other.
                double n1[3], n2[3], n12[3];
                findNormal(this->mesh.getPoint(a),this->mesh.getPoint(b),this->mesh.getPoint(c),n1);
                findNormal(this->mesh.getPoint(a),this->mesh.getPoint(b),this->mesh.getPoint(d),n2);
                cross(n1,n2,n12);
                if (dot(n1,n2) < 0.0 && length(n12) <= this->tolerance * length(n1) * length(n2))
                {
                    isSegment = false;
                }
            }
        }
        if (isSegment)
        {
            this->addSubSegment(a,b,noMarker);
        }
        else
        {
            facetEdges.push_back(key);
        }
    }

    // Order of edges in a hash table is not defined.
    std::sort(facetEdges.begin(),facetEdges.end());
    this->restoreDelaunay(facetEdges);

    for (size_t i=0;i<this->segments.size();i++)
    {
        this->queueSubSegment(int(i));
    }
    for (size_t i=0;i<this->faces.size();i++)
    {
        this->queueSubFace(int(i));
    }
}

bool RBoundaryRecovery::flipEdge(int a, int b, FlipMode mode, int &c, int &d)
{
    const uint64_t key = findEdgeKey(a,b);
    if (this->segmentMap.find(key) != this->segmentMap.end())
    {
        return false;
    }
    auto iter = this->edgeFaces.find(key);
    if (iter == this->edgeFaces.end() || iter->second.size() != 2)
    {
        return false;
    }
    const int faceID1 = iter->second[0];
    const int faceID2 = iter->second[1];
    const SubFace face1 = this->faces[size_t(faceID1)];
    const SubFace face2 = this->faces[size_t(faceID2)];

    // A recovered pair already agrees with the volume triangulation. Keep it,
    // including constrained pairs which are no longer Delaunay near sharp edges.
    if (mode == FlipDelaunay)
    {
        RTetraMesh::Face meshFace;
        if (this->mesh.findFace(face1.v[0],face1.v[1],face1.v[2],meshFace)
            && this->mesh.findFace(face2.v[0],face2.v[1],face2.v[2],meshFace))
        {
            return false;
        }
    }

    c = d = -1;
    for (int i=0;i<3;i++)
    {
        if (face1.v[i] != a && face1.v[i] != b)
        {
            c = face1.v[i];
        }
        if (face2.v[i] != a && face2.v[i] != b)
        {
            d = face2.v[i];
        }
    }
    if (c < 0 || d < 0 || c == d)
    {
        return false;
    }

    const double *pa = this->mesh.getPoint(a);
    const double *pb = this->mesh.getPoint(b);
    const double *pc = this->mesh.getPoint(c);
    const double *pd = this->mesh.getPoint(d);

    double n1[3], n2[3];
    findNormal(this->mesh.getPoint(face1.v[0]),this->mesh.getPoint(face1.v[1]),this->mesh.getPoint(face1.v[2]),n1);
    findNormal(this->mesh.getPoint(face2.v[0]),this->mesh.getPoint(face2.v[1]),this->mesh.getPoint(face2.v[2]),n2);
    const double n1Length = length(n1);
    if (!(n1Length > 0.0))
    {
        return false;
    }

    // The new edge has to cross the old one.
    double cd[3], ca[3], cb[3], cda[3], cdb[3];
    subtract(pd,pc,cd);
    subtract(pa,pc,ca);
    subtract(pb,pc,cb);
    cross(cd,ca,cda);
    cross(cd,cb,cdb);
    const double sideA = dot(cda,n1);
    const double sideB = dot(cdb,n1);
    const double cdLength = length(cd);
    if (!(sideA * sideB < 0.0)
        || std::fabs(sideA) <= 1.0e-10 * n1Length * cdLength * length(ca)
        || std::fabs(sideB) <= 1.0e-10 * n1Length * cdLength * length(cb))
    {
        return false;
    }

    double center[3], radius;
    if (!findCircumcircle(pa,pb,pc,center,radius))
    {
        return false;
    }
    const double dDistance = distance(pd,center);
    const double offset = (dDistance * dDistance - radius * radius) / (radius * radius);
    if (mode == FlipDelaunay)
    {
        if (!(offset < -flipTolerance))
        {
            return false;
        }
    }
    else
    {
        if (!(std::fabs(offset) < cocircularTolerance))
        {
            return false;
        }
    }

    // New subfaces keep orientation of the old ones.
    int t1[3] = {c,d,a};
    int t2[3] = {d,c,b};
    double nt[3];
    findNormal(this->mesh.getPoint(t1[0]),this->mesh.getPoint(t1[1]),this->mesh.getPoint(t1[2]),nt);
    if (dot(nt,n1) < 0.0)
    {
        std::swap(t1[0],t1[1]);
    }
    findNormal(this->mesh.getPoint(t2[0]),this->mesh.getPoint(t2[1]),this->mesh.getPoint(t2[2]),nt);
    if (dot(nt,n2) < 0.0)
    {
        std::swap(t2[0],t2[1]);
    }

    if (this->findSubFace(t1[0],t1[1],t1[2]) >= 0 || this->findSubFace(t2[0],t2[1],t2[2]) >= 0)
    {
        return false;
    }

    this->removeSubFace(faceID1);
    this->removeSubFace(faceID2);
    this->queueSubFace(this->addSubFace(t1[0],t1[1],t1[2],face1.marker));
    this->queueSubFace(this->addSubFace(t2[0],t2[1],t2[2],face2.marker));

    return true;
}

void RBoundaryRecovery::restoreDelaunay(std::vector<uint64_t> &edgeStack)
{
    // Number of flips is limited, as rounding errors could make them repeat.
    uint64_t nFlipsLeft = 100 * uint64_t(edgeStack.size()) + 10000;

    while (!edgeStack.empty() && nFlipsLeft > 0)
    {
        const uint64_t key = edgeStack.back();
        edgeStack.pop_back();

        const int a = int(key >> 32);
        const int b = int(key & 0xffffffffULL);
        int c, d;
        if (this->flipEdge(a,b,FlipDelaunay,c,d))
        {
            nFlipsLeft--;
            edgeStack.push_back(findEdgeKey(a,c));
            edgeStack.push_back(findEdgeKey(c,b));
            edgeStack.push_back(findEdgeKey(b,d));
            edgeStack.push_back(findEdgeKey(d,a));
        }
    }
    edgeStack.clear();
}

bool RBoundaryRecovery::matchSubFace(int faceID)
{
    const SubFace face = this->faces[size_t(faceID)];
    for (int i=0;i<3;i++)
    {
        const int a = face.v[i];
        const int b = face.v[(i+1)%3];
        const int c = face.v[(i+2)%3];
        const uint64_t key = findEdgeKey(a,b);
        if (this->segmentMap.find(key) != this->segmentMap.end())
        {
            continue;
        }
        auto iter = this->edgeFaces.find(key);
        if (iter == this->edgeFaces.end() || iter->second.size() != 2)
        {
            continue;
        }
        const SubFace &rOther = this->faces[size_t(iter->second[0] == faceID ? iter->second[1] : iter->second[0])];
        int d = -1;
        for (int j=0;j<3;j++)
        {
            if (rOther.v[j] != a && rOther.v[j] != b)
            {
                d = rOther.v[j];
            }
        }
        if (d < 0 || d == c)
        {
            continue;
        }
        RTetraMesh::Face meshFace;
        if (!this->mesh.findFace(c,d,a,meshFace) || !this->mesh.findFace(c,d,b,meshFace))
        {
            continue;
        }
        int c2, d2;
        if (this->flipEdge(a,b,FlipCocircular,c2,d2))
        {
            return true;
        }
    }
    return false;
}

bool RBoundaryRecovery::matchFacet(int faceID)
{
    // Reuse the volume mesh's triangulation instead of repeatedly subdividing
    // a planar facet to force a different choice of diagonals. Nothing changes
    // until the entire replacement has been validated.
    const int pointCount = this->mesh.getNPoints();
    if (this->facetMatchPoints[size_t(faceID)] == pointCount) return false;
    const SubFace seed = this->faces[size_t(faceID)];
    const double *origin = this->mesh.getPoint(seed.v[0]);
    double normal[3];
    findNormal(origin,this->mesh.getPoint(seed.v[1]),this->mesh.getPoint(seed.v[2]),normal);
    const double normalLength = length(normal);
    if (!(normalLength > 0.0)) return false;
    for (double &component : normal) component /= normalLength;

    const unsigned int stamp = ++this->facetStamp;
    std::vector<int> patch{faceID};
    this->facetMarks[size_t(faceID)] = stamp;
    std::set<int> vertices;
    std::map<uint64_t,int> oldEdges;
    std::set<uint64_t> constraints;
    double oldArea = 0.0;
    for (size_t i=0;i<patch.size();i++)
    {
        // Bound a fallback independently of the Steiner-node budget.
        if (patch.size() > 4096) return false;
        const int id = patch[i];
        this->facetMatchPoints[size_t(id)] = pointCount;
        const SubFace &face = this->faces[size_t(id)];
        double n[3];
        findNormal(this->mesh.getPoint(face.v[0]),this->mesh.getPoint(face.v[1]),this->mesh.getPoint(face.v[2]),n);
        if (!(dot(n,normal) > 0.0)) return false;
        oldArea += dot(n,normal);
        for (int j=0;j<3;j++)
        {
            vertices.insert(face.v[j]);
            const uint64_t key = findEdgeKey(face.v[j],face.v[(j+1)%3]);
            oldEdges[key]++;
            if (this->segmentMap.count(key))
            {
                constraints.insert(key);
                continue;
            }
            auto adjacent = this->edgeFaces.find(key);
            if (adjacent == this->edgeFaces.end() || adjacent->second.size() != 2) return false;
            for (int neighbor : adjacent->second)
            {
                if (this->facetMarks[size_t(neighbor)] == stamp) continue;
                const SubFace &other = this->faces[size_t(neighbor)];
                if (other.marker != seed.marker) return false;
                // Do not accumulate approximate coplanarity across a curved patch.
                for (int v : other.v)
                {
                    double offset[3];subtract(this->mesh.getPoint(v),origin,offset);
                    if (std::fabs(dot(offset,normal)) > 64*std::numeric_limits<double>::epsilon()
                        * std::max(1.0,length(offset))) return false;
                }
                this->facetMarks[size_t(neighbor)] = stamp;
                patch.push_back(neighbor);
            }
        }
    }

    auto inside = [&](const double *p) {
        for (int id : patch)
        {
            const SubFace &face = this->faces[size_t(id)];
            bool contains = true;
            for (int j=0;j<3;j++)
            {
                const double *a = this->mesh.getPoint(face.v[j]);
                double edge[3], offset[3], product[3];
                subtract(this->mesh.getPoint(face.v[(j+1)%3]),a,edge);
                subtract(p,a,offset);cross(edge,offset,product);
                if (dot(product,normal) < -64*std::numeric_limits<double>::epsilon()*length(edge)*length(offset))
                {
                    contains = false;
                    break;
                }
            }
            if (contains) return true;
        }
        return false;
    };

    std::set<std::array<int,3>> seen;
    std::vector<std::array<int,3>> replacement;
    std::map<uint64_t,int> newEdges;
    std::set<int> newVertices;
    double newArea = 0.0;
    for (int vertex : vertices)
    {
        std::vector<int> star;
        this->mesh.findStar(vertex,star);
        for (int tetraID : star)
        {
            const auto &tetra = this->mesh.getTetra(tetraID);
            for (int j=0;j<4;j++)
            {
                const int *indices = RTetraMesh::faceVertices[j];
                std::array<int,3> triangle{tetra.v[indices[0]],tetra.v[indices[1]],tetra.v[indices[2]]};
                if (!vertices.count(triangle[0]) || !vertices.count(triangle[1]) || !vertices.count(triangle[2])) continue;
                auto key = triangle;std::sort(key.begin(),key.end());
                if (!seen.insert(key).second) continue;
                const double *a = this->mesh.getPoint(triangle[0]);
                const double *b = this->mesh.getPoint(triangle[1]);
                const double *c = this->mesh.getPoint(triangle[2]);
                const double center[3] = {a[0]+((b[0]-a[0])+(c[0]-a[0]))/3,
                                          a[1]+((b[1]-a[1])+(c[1]-a[1]))/3,
                                          a[2]+((b[2]-a[2])+(c[2]-a[2]))/3};
                if (!inside(center)) continue;
                double n[3];findNormal(a,b,c,n);
                const double area = dot(n,normal);
                if (area == 0.0) return false;
                if (area < 0.0) std::swap(triangle[1],triangle[2]);
                newArea += std::fabs(area);
                for (int k=0;k<3;k++)
                {
                    newVertices.insert(triangle[k]);
                    newEdges[findEdgeKey(triangle[k],triangle[(k+1)%3])]++;
                }
                replacement.push_back(triangle);
            }
        }
    }
    if (replacement.empty() || newVertices != vertices
        || std::fabs(newArea-oldArea) > 1e-10*oldArea) return false;
    for (const auto &edge : oldEdges)
    {
        if (edge.second > 2 || (edge.second == 1 && newEdges[edge.first] != 1)) return false;
    }
    for (const auto &edge : newEdges)
    {
        if (edge.second > 2 || (edge.second == 1 && oldEdges[edge.first] != 1)) return false;
    }
    for (uint64_t edge : constraints) if (!newEdges.count(edge)) return false;

    for (int id : patch) this->removeSubFace(id);
    for (const auto &triangle : replacement)
    {
        this->queueSubFace(this->addSubFace(triangle[0],triangle[1],triangle[2],seed.marker));
    }
    return true;
}

RBoundaryRecovery::FacetLocation RBoundaryRecovery::locateInFacet(int faceID, const double *p) const
{
    FacetLocation location;

    int currentID = faceID;
    uint64_t previousKey = std::numeric_limits<uint64_t>::max();

    for (int step=0;step<1000000;step++)
    {
        location.visitedFaceIDs.push_back(currentID);

        const SubFace &rFace = this->faces[size_t(currentID)];
        const double *pv[3] = {this->mesh.getPoint(rFace.v[0]),this->mesh.getPoint(rFace.v[1]),this->mesh.getPoint(rFace.v[2])};

        double n[3];
        findNormal(pv[0],pv[1],pv[2],n);
        const double nLength = length(n);
        if (!(nLength > 0.0))
        {
            return location;
        }

        // Signed distances of the point from edges, positive inside of the subface.
        double distances[3];
        double scale = 0.0;
        for (int i=0;i<3;i++)
        {
            double edge[3], toPoint[3], product[3];
            subtract(pv[(i+1)%3],pv[i],edge);
            subtract(p,pv[i],toPoint);
            cross(edge,toPoint,product);
            const double edgeLength = length(edge);
            scale = std::max(scale,edgeLength);
            distances[i] = dot(product,n) / (nLength * edgeLength);
        }
        const double toleranceDistance = onEdgeTolerance * scale;

        int worst = -1;
        double worstDistance = -toleranceDistance;
        for (int i=0;i<3;i++)
        {
            if (findEdgeKey(rFace.v[i],rFace.v[(i+1)%3]) == previousKey)
            {
                continue;
            }
            if (distances[i] < worstDistance)
            {
                worstDistance = distances[i];
                worst = i;
            }
        }

        if (worst < 0)
        {
            location.faceID = currentID;
            for (int i=0;i<3;i++)
            {
                if (std::fabs(distances[i]) <= toleranceDistance)
                {
                    location.onEdge = true;
                    location.edge[0] = rFace.v[i];
                    location.edge[1] = rFace.v[(i+1)%3];
                    location.segmentID = this->findSubSegment(location.edge[0],location.edge[1]);
                    break;
                }
            }
            return location;
        }

        const int a = rFace.v[worst];
        const int b = rFace.v[(worst+1)%3];
        const int segmentID = this->findSubSegment(a,b);
        if (segmentID >= 0)
        {
            location.segmentID = segmentID;
            return location;
        }
        const uint64_t key = findEdgeKey(a,b);
        auto iter = this->edgeFaces.find(key);
        if (iter == this->edgeFaces.end() || iter->second.size() != 2)
        {
            return location;
        }
        previousKey = key;
        currentID = iter->second[0] == currentID ? iter->second[1] : iter->second[0];
    }
    return location;
}

void RBoundaryRecovery::splitEdgeFaces(int a, int b, int vertexID, std::vector<uint64_t> &edgeStack)
{
    auto iter = this->edgeFaces.find(findEdgeKey(a,b));
    if (iter == this->edgeFaces.end())
    {
        return;
    }
    const std::vector<int> faceIDs = iter->second;
    for (int faceID : faceIDs)
    {
        const SubFace face = this->faces[size_t(faceID)];
        int first[3], second[3];
        int c = -1;
        for (int i=0;i<3;i++)
        {
            first[i] = face.v[i] == b ? vertexID : face.v[i];
            second[i] = face.v[i] == a ? vertexID : face.v[i];
            if (face.v[i] != a && face.v[i] != b)
            {
                c = face.v[i];
            }
        }
        this->removeSubFace(faceID);
        this->queueSubFace(this->addSubFace(first[0],first[1],first[2],face.marker));
        this->queueSubFace(this->addSubFace(second[0],second[1],second[2],face.marker));
        if (c >= 0)
        {
            edgeStack.push_back(findEdgeKey(a,c));
            edgeStack.push_back(findEdgeKey(b,c));
        }
    }
}

bool RBoundaryRecovery::splitSubSegment(int segmentID, double minRadius)
{
    this->lastSplitFailure = "Constrained edge could not be split while preserving the boundary";
    const SubSegment segment = this->segments[size_t(segmentID)];
    const int a = segment.v[0];
    const int b = segment.v[1];
    const double *pa = this->mesh.getPoint(a);
    const double *pb = this->mesh.getPoint(b);

    const double segmentLength = distance(pa,pb);
    if (!(segmentLength > this->minLength))
    {
        this->lastSplitFailure = "Constrained edge length is at or below the minimum split length";
        return false;
    }

    // Subsegments which meet in an input vertex are split in the same distance from it, a power of two,
    // so that their parts do not cause each other to be split again and again.
    const bool aIsInput = this->vertexTypes[size_t(a)] == VertexInput;
    const bool bIsInput = this->vertexTypes[size_t(b)] == VertexInput;
    double ratio = 0.5;
    if (aIsInput != bIsInput)
    {
        ratio = std::exp2(std::round(std::log2(0.5 * segmentLength))) / segmentLength;
        if (!(ratio > 0.2 && ratio < 0.8))
        {
            ratio = 0.5;
        }
        if (bIsInput)
        {
            ratio = 1.0 - ratio;
        }
    }
    const double p[3] = {pa[0] + ratio * (pb[0] - pa[0]),
                         pa[1] + ratio * (pb[1] - pa[1]),
                         pa[2] + ratio * (pb[2] - pa[2])};

    const double sizeA = this->vertexSizes[size_t(a)];
    const double sizeB = this->vertexSizes[size_t(b)];
    const uint64_t key = findEdgeKey(a,b);

    if (minRadius > 0.0)
    {
        // Subsegments which meet in an input vertex are split in the same distance from it. The new vertex
        // will therefore get a neighbor on each of them, which must not be too close either.
        std::vector<int> star;
        for (int end=0;end<2;end++)
        {
            const int vertexID = end == 0 ? a : b;
            const int otherID = end == 0 ? b : a;
            if (this->vertexTypes[size_t(vertexID)] != VertexInput)
            {
                continue;
            }
            const double *pv = this->mesh.getPoint(vertexID);
            double toPoint[3];
            subtract(p,pv,toPoint);
            const double pointDistance = length(toPoint);
            this->mesh.findStar(vertexID,star);
            for (int tetraID : star)
            {
                const RTetraMesh::Tetra &rTetra = this->mesh.getTetra(tetraID);
                for (int i=0;i<4;i++)
                {
                    const int neighborID = rTetra.v[i];
                    if (neighborID == RTetraMesh::ghostVertex || neighborID == vertexID || neighborID == otherID
                        || !this->boundaryVertices[size_t(neighborID)] || this->findSubSegment(vertexID,neighborID) < 0)
                    {
                        continue;
                    }
                    double toNeighbor[3];
                    subtract(this->mesh.getPoint(neighborID),pv,toNeighbor);
                    const double cosAngle = dot(toPoint,toNeighbor) / (pointDistance * length(toNeighbor));
                    if (cosAngle > 0.0 && pointDistance * std::sqrt(2.0 * std::max(0.0,1.0 - cosAngle)) < minRadius)
                    {
                        return false;
                    }
                }
            }
        }
    }

    const int locatedID = this->mesh.locate(p,this->mesh.getVertexTetra(a));
    const int duplicateID = this->mesh.findDuplicateVertex(locatedID,p);
    if (duplicateID >= 0)
    {
        if (duplicateID == a || duplicateID == b)
        {
            this->lastSplitFailure = "Edge split point rounds to an existing endpoint";
            return false;
        }
        // There is a vertex lying on the subsegment already, which splits it.
        this->removeSubSegment(segmentID);
        this->queueSubSegment(this->addSubSegment(a,duplicateID,segment.marker));
        this->queueSubSegment(this->addSubSegment(duplicateID,b,segment.marker));
        std::vector<uint64_t> edgeStack;
        this->splitEdgeFaces(a,b,duplicateID,edgeStack);
        this->restoreDelaunay(edgeStack);
        return true;
    }

    std::vector<int> incidentFaceIDs;
    {
        auto iter = this->edgeFaces.find(key);
        if (iter != this->edgeFaces.end())
        {
            incidentFaceIDs = iter->second;
        }
    }

    // Subfaces which are replaced: at first all those which circumscribed circles contain the new vertex,
    // then only those at the subsegment, which can always be replaced.
    std::vector<int> removedFaceIDs;
    bool cavityFound = false;
    for (int attempt=0;attempt<2 && !cavityFound;attempt++)
    {
        removedFaceIDs.clear();
        bool valid = true;
        for (int faceID : incidentFaceIDs)
        {
            std::vector<int> facetFaceIDs;
            if (attempt == 0)
            {
                this->findFacetCavity(p,std::vector<int>{faceID},key,facetFaceIDs);
                valid = valid && !facetFaceIDs.empty();
            }
            else
            {
                facetFaceIDs.push_back(faceID);
            }
            removedFaceIDs.insert(removedFaceIDs.end(),facetFaceIDs.begin(),facetFaceIDs.end());
        }
        cavityFound = valid && this->findBoundaryCavity(p,removedFaceIDs,a,b,locatedID,minRadius == 0.0);
    }
    if (!cavityFound)
    {
        this->lastSplitFailure = "No boundary-preserving cavity could be constructed for the constrained edge";
        return false;
    }
    if (minRadius > 0.0 && !this->insertionAllowed(p,VertexSegment,minRadius,{a,b}))
    {
        return false;
    }

    const std::vector<LoopEdge> loopEdges = this->findLoopEdges(removedFaceIDs,key);

    const int vertexID = this->insertCavityPoint(p,VertexSegment,minRadius,{a,b});
    if (sizeA > 0.0 && sizeB > 0.0)
    {
        this->vertexSizes[size_t(vertexID)] = sizeA + ratio * (sizeB - sizeA);
    }

    this->removeSubSegment(segmentID);
    this->queueSubSegment(this->addSubSegment(a,vertexID,segment.marker));
    this->queueSubSegment(this->addSubSegment(vertexID,b,segment.marker));

    for (int faceID : removedFaceIDs)
    {
        this->removeSubFace(faceID);
    }
    std::vector<uint64_t> edgeStack;
    for (const LoopEdge &rEdge : loopEdges)
    {
        this->queueSubFace(this->addSubFace(rEdge.a,rEdge.b,vertexID,rEdge.marker));
        edgeStack.push_back(findEdgeKey(rEdge.a,rEdge.b));
    }
    // Splitting must preserve the facet's constrained Delaunay triangulation.
    // Otherwise a stale circumcenter can repeatedly split a distant subsegment.
    this->restoreDelaunay(edgeStack);
    return true;
}

bool RBoundaryRecovery::splitSubFace(int faceID, double minRadius)
{
    this->lastSplitFailure = "No valid boundary cavity could be constructed";
    if (this->splitSubFaceAt(faceID,minRadius,false)) return true;
    // Recovery may use an interior point when the circumcenter hits a protected
    // edge, an existing vertex, or an invalid cavity. Refinement keeps its sizing
    // and encroachment safeguards and does not take this recovery-only fallback.
    return minRadius == 0.0 && this->splitSubFaceAt(faceID,minRadius,true);
}

bool RBoundaryRecovery::splitSubFaceAt(int faceID, double minRadius, bool centroid)
{
    auto fail = [&](const char *reason) { this->lastSplitFailure = reason;return false; };
    const SubFace face = this->faces[size_t(faceID)];

    double center[3], radius;
    if (!findCircumcircle(this->mesh.getPoint(face.v[0]),this->mesh.getPoint(face.v[1]),this->mesh.getPoint(face.v[2]),center,radius))
    {
        return fail("Triangle has a degenerate or non-finite circumcircle");
    }
    if (!(radius > this->minLength))
    {
        return fail("Triangle circumradius is at or below the minimum split length");
    }

    if (centroid)
    {
        for (int i=0;i<3;i++)
        {
            const double a = this->mesh.getPoint(face.v[0])[i];
            center[i] = a+((this->mesh.getPoint(face.v[1])[i]-a)+(this->mesh.getPoint(face.v[2])[i]-a))/3.0;
        }
    }
    const FacetLocation location = this->locateInFacet(faceID,center);
    if (location.segmentID >= 0)
    {
        // The center lies behind a subsegment or on it.
        return this->splitSubSegment(location.segmentID,minRadius);
    }
    if (location.faceID < 0)
    {
        return fail("Split point could not be located in the surface facet");
    }

    // During quality refinement protect subsegments' diametral spheres. Recovery
    // only needs to preserve the segments themselves (checked by findLostSubSegment).
    // Enforcing empty spheres here can cause endless splitting at facet corners.
    for (int visitedID : location.visitedFaceIDs)
    {
        const SubFace &rVisited = this->faces[size_t(visitedID)];
        for (int i=0;i<3;i++)
        {
            const int a = rVisited.v[i];
            const int b = rVisited.v[(i+1)%3];
            const int segmentID = this->findSubSegment(a,b);
            if (minRadius > 0.0 && segmentID >= 0 && isInDiametralSphere(center,this->mesh.getPoint(a),this->mesh.getPoint(b)))
            {
                return this->splitSubSegment(segmentID,minRadius);
            }
        }
    }

    // Subfaces of a facet do not have to lie exactly in one plane.
    const SubFace target = this->faces[size_t(location.faceID)];
    double p[3] = {center[0],center[1],center[2]};
    {
        double n[3], toPoint[3];
        findNormal(this->mesh.getPoint(target.v[0]),this->mesh.getPoint(target.v[1]),this->mesh.getPoint(target.v[2]),n);
        subtract(center,this->mesh.getPoint(target.v[0]),toPoint);
        const double n2 = dot(n,n);
        if (n2 > 0.0)
        {
            const double offset = dot(toPoint,n) / n2;
            for (int i=0;i<3;i++)
            {
                p[i] = center[i] - offset * n[i];
            }
        }
    }

    if (location.onEdge)
    {
        const double *a=this->mesh.getPoint(location.edge[0]), *b=this->mesh.getPoint(location.edge[1]);
        double edge[3],offset[3];subtract(b,a,edge);subtract(center,a,offset);
        const double ratio=dot(offset,edge)/dot(edge,edge);
        if (!(ratio>0.0 && ratio<1.0)) return fail("Split point is not strictly inside the facet edge");
        for (int i=0;i<3;i++) p[i]=a[i]+ratio*edge[i];
    }

    // Subfaces which the new vertex lies in.
    std::vector<int> startFaceIDs;
    if (location.onEdge)
    {
        auto iter = this->edgeFaces.find(findEdgeKey(location.edge[0],location.edge[1]));
        if (iter == this->edgeFaces.end() || iter->second.size() != 2)
        {
            return fail("Facet edge does not have exactly two incident surface triangles");
        }
        startFaceIDs = iter->second;
    }
    else
    {
        startFaceIDs.push_back(location.faceID);
    }

    const int locatedID = this->mesh.locate(p,this->mesh.getVertexTetra(target.v[0]));
    const int duplicateID = this->mesh.findDuplicateVertex(locatedID,p);
    if (duplicateID >= 0)
    {
        // There is a vertex in the subface already, it becomes a part of the facet.
        for (int startID : startFaceIDs)
        {
            const SubFace &rStart = this->faces[size_t(startID)];
            if (duplicateID == rStart.v[0] || duplicateID == rStart.v[1] || duplicateID == rStart.v[2])
            {
                return fail("Split point rounds to an existing triangle vertex");
            }
        }
        std::vector<uint64_t> edgeStack;
        if (location.onEdge)
        {
            this->splitEdgeFaces(location.edge[0],location.edge[1],duplicateID,edgeStack);
        }
        else
        {
            this->removeSubFace(location.faceID);
            for (int i=0;i<3;i++)
            {
                this->queueSubFace(this->addSubFace(target.v[i],target.v[(i+1)%3],duplicateID,target.marker));
                edgeStack.push_back(findEdgeKey(target.v[i],target.v[(i+1)%3]));
            }
        }
        this->restoreDelaunay(edgeStack);
        return true;
    }

    // Subfaces which are replaced: at first all those which circumscribed circles contain the new vertex,
    // then only those which the vertex lies in, which can always be replaced.
    const uint64_t noKey = std::numeric_limits<uint64_t>::max();
    std::vector<int> removedFaceIDs;
    std::vector<LoopEdge> loopEdges;
    bool cavityFound = false;
    for (int attempt=0;attempt<2 && !cavityFound;attempt++)
    {
        if (attempt == 0)
        {
            this->findFacetCavity(p,startFaceIDs,noKey,removedFaceIDs);
            if (removedFaceIDs.empty())
            {
                continue;
            }
        }
        else
        {
            removedFaceIDs = startFaceIDs;
        }
        loopEdges = this->findLoopEdges(removedFaceIDs,noKey);
        for (const LoopEdge &rEdge : loopEdges)
        {
            const int segmentID = this->findSubSegment(rEdge.a,rEdge.b);
            if (minRadius > 0.0 && segmentID >= 0 && isInDiametralSphere(p,this->mesh.getPoint(rEdge.a),this->mesh.getPoint(rEdge.b)))
            {
                return this->splitSubSegment(segmentID,minRadius);
            }
        }
        cavityFound = this->findBoundaryCavity(p,removedFaceIDs,-1,-1,locatedID,minRadius == 0.0);
    }
    if (!cavityFound)
    {
        return fail("No boundary-preserving cavity could be constructed at the split point");
    }

    // A subsegment which would be removed is split instead.
    const int lostSegmentID = this->findLostSubSegment(-1);
    if (lostSegmentID >= 0)
    {
        return this->splitSubSegment(lostSegmentID,minRadius);
    }
    if (minRadius > 0.0 && !this->insertionAllowed(p,VertexFacet,minRadius,{target.v[0],target.v[1],target.v[2]}))
    {
        return fail("Split would violate the refinement spacing safeguard");
    }

    const int vertexID = this->insertCavityPoint(p,VertexFacet,minRadius,{target.v[0],target.v[1],target.v[2]});

    for (int removedID : removedFaceIDs)
    {
        this->removeSubFace(removedID);
    }
    std::vector<uint64_t> edgeStack;
    for (const LoopEdge &rEdge : loopEdges)
    {
        this->queueSubFace(this->addSubFace(rEdge.a,rEdge.b,vertexID,rEdge.marker));
        edgeStack.push_back(findEdgeKey(rEdge.a,rEdge.b));
    }
    // Splitting must preserve the facet's constrained Delaunay triangulation.
    // Otherwise a stale circumcenter can repeatedly split a distant subsegment.
    this->restoreDelaunay(edgeStack);
    return true;
}

void RBoundaryRecovery::recoverSubFace(int faceID)
{
    if (this->matchSubFace(faceID) || this->matchFacet(faceID))
    {
        return;
    }

    // The subface is checked again, as it may still be missing after a vertex is added.
    this->queueSubFace(faceID);

    if (!this->splitSubFace(faceID))
    {
        const SubFace &rFace = this->faces[size_t(faceID)];
        const double *pa = this->mesh.getPoint(rFace.v[0]);
        const double *pb = this->mesh.getPoint(rFace.v[1]);
        const double *pc = this->mesh.getPoint(rFace.v[2]);
        double center[3], radius = std::numeric_limits<double>::quiet_NaN();
        findCircumcircle(pa,pb,pc,center,radius);
        throw RError(RError::Type::Application,R_ERROR_REF,
                     "Failed to recover boundary triangle %d (nodes %d, %d, %d, surface marker %d): %s. "
                     "Circumradius = %.17g, minimum split length = %.17g. "
                     "Coordinates: (%.17g %.17g %.17g), (%.17g %.17g %.17g), (%.17g %.17g %.17g). "
                     "Boundary recovery did not converge; this does not establish a surface intersection.",
                     faceID,rFace.v[0],rFace.v[1],rFace.v[2],rFace.marker,this->lastSplitFailure,radius,this->minLength,
                     pa[0],pa[1],pa[2],pb[0],pb[1],pb[2],pc[0],pc[1],pc[2]);
    }
}

void RBoundaryRecovery::recover()
{
    while (true)
    {
        // Subsegments go first, as subfaces can be recovered only between them.
        if (!this->segmentQueue.empty())
        {
            const int segmentID = this->segmentQueue.back();
            this->segmentQueue.pop_back();
            SubSegment &rSegment = this->segments[size_t(segmentID)];
            rSegment.queued = false;
            if (rSegment.alive && this->mesh.findEdge(rSegment.v[0],rSegment.v[1]) == RTetraMesh::noTetra)
            {
                const int a = rSegment.v[0];
                const int b = rSegment.v[1];
                if (!this->splitSubSegment(segmentID))
                {
                    const double *pa = this->mesh.getPoint(a);
                    const double *pb = this->mesh.getPoint(b);
                    throw RError(RError::Type::Application,R_ERROR_REF,
                                 "Failed to recover boundary edge %d (nodes %d, %d, marker %d): %s. "
                                 "Edge length = %.17g, minimum split length = %.17g. "
                                 "Coordinates: (%.17g %.17g %.17g), (%.17g %.17g %.17g). "
                                 "Boundary recovery did not converge; this does not establish a surface intersection.",
                                 segmentID,a,b,rSegment.marker,this->lastSplitFailure,distance(pa,pb),this->minLength,
                                 pa[0],pa[1],pa[2],pb[0],pb[1],pb[2]);
                }
            }
            continue;
        }
        if (!this->faceQueue.empty())
        {
            const int faceID = this->faceQueue.back();
            this->faceQueue.pop_back();
            SubFace &rFace = this->faces[size_t(faceID)];
            rFace.queued = false;
            RTetraMesh::Face meshFace;
            if (rFace.alive && !this->mesh.findFace(rFace.v[0],rFace.v[1],rFace.v[2],meshFace))
            {
                this->recoverSubFace(faceID);
            }
            continue;
        }
        break;
    }
}

const std::vector<RBoundaryRecovery::SubSegment> &RBoundaryRecovery::getSubSegments() const
{
    return this->segments;
}

const std::vector<RBoundaryRecovery::SubFace> &RBoundaryRecovery::getSubFaces() const
{
    return this->faces;
}

int RBoundaryRecovery::findSubSegment(int a, int b) const
{
    auto iter = this->segmentMap.find(findEdgeKey(a,b));
    return iter == this->segmentMap.end() ? -1 : iter->second;
}

int RBoundaryRecovery::findSubFace(int a, int b, int c) const
{
    FaceKey key = {{a,b,c}};
    std::sort(key.v,key.v+3);
    auto iter = this->faceMap.find(key);
    return iter == this->faceMap.end() ? -1 : iter->second;
}

int RBoundaryRecovery::findSubFace(const RTetraMesh::Face &face) const
{
    const RTetraMesh::Tetra &rTetra = this->mesh.getTetra(face.tetra);
    const int *pFace = RTetraMesh::faceVertices[face.face];
    const int a = rTetra.v[pFace[0]];
    const int b = rTetra.v[pFace[1]];
    const int c = rTetra.v[pFace[2]];
    if (a == RTetraMesh::ghostVertex || b == RTetraMesh::ghostVertex || c == RTetraMesh::ghostVertex)
    {
        return -1;
    }
    if (!this->boundaryVertices[size_t(a)] || !this->boundaryVertices[size_t(b)] || !this->boundaryVertices[size_t(c)])
    {
        return -1;
    }
    return this->findSubFace(a,b,c);
}

double RBoundaryRecovery::findSubSegmentLength(int segmentID) const
{
    const SubSegment &rSegment = this->segments[size_t(segmentID)];
    return distance(this->mesh.getPoint(rSegment.v[0]),this->mesh.getPoint(rSegment.v[1]));
}

double RBoundaryRecovery::findSubFaceRadius(int faceID) const
{
    const SubFace &rFace = this->faces[size_t(faceID)];
    double center[3], radius;
    if (!findCircumcircle(this->mesh.getPoint(rFace.v[0]),this->mesh.getPoint(rFace.v[1]),this->mesh.getPoint(rFace.v[2]),center,radius))
    {
        return 0.0;
    }
    return radius;
}

RBoundaryRecovery::Encroachment RBoundaryRecovery::findEncroachment(const double *p) const
{
    Encroachment encroachment;

    // Subsegment which the point lies too close to.
    for (int tetraID : this->mesh.getCavityTetras())
    {
        const RTetraMesh::Tetra &rTetra = this->mesh.getTetra(tetraID);
        for (int i=0;i<6;i++)
        {
            const int a = rTetra.v[RTetraMesh::edgeVertices[i][0]];
            const int b = rTetra.v[RTetraMesh::edgeVertices[i][1]];
            if (a == RTetraMesh::ghostVertex || b == RTetraMesh::ghostVertex)
            {
                continue;
            }
            if (!this->boundaryVertices[size_t(a)] || !this->boundaryVertices[size_t(b)])
            {
                continue;
            }
            const int segmentID = this->findSubSegment(a,b);
            if (segmentID >= 0 && isInDiametralSphere(p,this->mesh.getPoint(a),this->mesh.getPoint(b)))
            {
                encroachment.type = EncroachmentSegment;
                encroachment.id = segmentID;
                return encroachment;
            }
        }
    }

    // Subsegment which would be removed.
    const int lostSegmentID = this->findLostSubSegment(-1);
    if (lostSegmentID >= 0)
    {
        encroachment.type = EncroachmentSegment;
        encroachment.id = lostSegmentID;
        return encroachment;
    }

    // Subface which the point lies too close to. Subfaces can be found only on boundary of the cavity.
    for (const RTetraMesh::Face &rMeshFace : this->mesh.getCavityFaces())
    {
        const int faceID = this->findSubFace(rMeshFace);
        if (faceID < 0)
        {
            continue;
        }
        const SubFace &rFace = this->faces[size_t(faceID)];
        double center[3], radius;
        if (findCircumcircle(this->mesh.getPoint(rFace.v[0]),this->mesh.getPoint(rFace.v[1]),this->mesh.getPoint(rFace.v[2]),center,radius)
            && distance(p,center) < radius)
        {
            encroachment.type = EncroachmentFace;
            encroachment.id = faceID;
            return encroachment;
        }
    }

    return encroachment;
}

int RBoundaryRecovery::classifyRegions()
{
    // Region of tetrahedra which were not visited yet.
    const int unknownRegion = -2;

    const int nSlots = this->mesh.getNTetraSlots();
    std::vector<int> stack;

    for (int i=0;i<nSlots;i++)
    {
        if (!this->mesh.isAlive(i))
        {
            continue;
        }
        if (this->mesh.isGhost(i))
        {
            this->mesh.setRegion(i,RTetraMesh::noRegion);
            stack.push_back(i);
        }
        else
        {
            this->mesh.setRegion(i,unknownRegion);
        }
    }

    auto flood = [&](int region)
    {
        while (!stack.empty())
        {
            const int tetraID = stack.back();
            stack.pop_back();
            for (int i=0;i<4;i++)
            {
                const int neighborID = this->mesh.getTetra(tetraID).n[i];
                if (this->mesh.getTetra(neighborID).region != unknownRegion)
                {
                    continue;
                }
                if (this->findSubFace(RTetraMesh::Face{tetraID,i}) >= 0)
                {
                    continue;
                }
                this->mesh.setRegion(neighborID,region);
                stack.push_back(neighborID);
            }
        }
    };

    flood(RTetraMesh::noRegion);

    int nRegions = 0;
    for (int i=0;i<nSlots;i++)
    {
        if (!this->mesh.isAlive(i) || this->mesh.getTetra(i).region != unknownRegion)
        {
            continue;
        }
        this->mesh.setRegion(i,nRegions);
        stack.push_back(i);
        flood(nRegions);
        nRegions++;
    }
    return nRegions;
}

void RBoundaryRecovery::setCollectCreatedTetras(bool collectCreatedTetras)
{
    this->collectCreatedTetras = collectCreatedTetras;
    if (!collectCreatedTetras)
    {
        this->createdTetras.clear();
    }
}

std::vector<int> &RBoundaryRecovery::getCreatedTetras()
{
    return this->createdTetras;
}

bool RBoundaryRecovery::insertionAllowed(const double *p, VertexType type, double minRadius, const std::vector<int> &parents) const
{
    std::vector<int> ancestry;
    for (int parent : parents)
        ancestry.insert(ancestry.end(),this->vertexParents[size_t(parent)].begin(),this->vertexParents[size_t(parent)].end());
    std::sort(ancestry.begin(),ancestry.end());
    for (int id : this->mesh.getCavityTetras()) for (int vertex : this->mesh.getTetra(id).v)
    {
        if (vertex==RTetraMesh::ghostVertex) continue;
        bool adjacent=type==VertexVolume || this->vertexTypes[size_t(vertex)]==VertexVolume;
        for (int root : this->vertexParents[size_t(vertex)])
            adjacent=adjacent || std::binary_search(ancestry.begin(),ancestry.end(),root);
        const double limit=adjacent?std::max(minRadius,this->insertionRadii[size_t(vertex)]):minRadius;
        if (distance(p,this->mesh.getPoint(vertex))<limit*(1.0-1e-10)) return false;
    }
    return true;
}

RBoundaryRecovery::RBoundaryRecovery(RTetraMesh &mesh, const RBoundaryRecovery &source)
    : RBoundaryRecovery(mesh,source.tolerance)
{
    this->copyState(source);
}

void RBoundaryRecovery::copyState(const RBoundaryRecovery &source)
{
    this->tolerance=source.tolerance;
    this->minLength=source.minLength;
    this->maxSteinerVertices=source.maxSteinerVertices;
    this->segments=source.segments;
    this->faces=source.faces;
    this->segmentMap=source.segmentMap;
    this->faceMap=source.faceMap;
    this->edgeFaces=source.edgeFaces;
    this->segmentQueue=source.segmentQueue;
    this->faceQueue=source.faceQueue;
    this->vertexTypes=source.vertexTypes;
    this->boundaryVertices=source.boundaryVertices;
    this->vertexSizes=source.vertexSizes;
    this->insertionRadii=source.insertionRadii;
    this->vertexParents=source.vertexParents;
    this->replacedMarks=source.replacedMarks;
    this->replacedStamp=source.replacedStamp;
    this->facetMarks=source.facetMarks;
    this->facetMatchPoints=source.facetMatchPoints;
    this->lastSplitFailure=source.lastSplitFailure;
    this->facetStamp=source.facetStamp;
    this->collectCreatedTetras=source.collectCreatedTetras;
    this->createdTetras=source.createdTetras;
    std::copy(source.nSteinerVertices,source.nSteinerVertices+nVertexTypes,this->nSteinerVertices);
    this->blockedFunction=[this](int tetra,int face) { return this->isFaceBlocked(tetra,face); };
}
