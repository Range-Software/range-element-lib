#ifndef REL_BOUNDARY_RECOVERY_H
#define REL_BOUNDARY_RECOVERY_H

#include <cstdint>
#include <functional>
#include <unordered_map>
#include <vector>

#include "rel_tetra_mesh.h"

//! Boundary of a meshed domain and its recovery in a tetrahedralization.
//!
//! The boundary is made of subfaces (triangles) and subsegments (edges which must not be removed:
//! lines, borders of surfaces and edges where surfaces meet or bend). Neighboring subfaces which
//! lie in one plane, belong to the same surface and are not separated by a subsegment form a facet.
//! Triangulation of a facet is free to change, only its shape is kept.
//!
//! A subsegment or a subface which is missing in the tetrahedralization is split by a new vertex,
//! until all of them are present. Subfaces which are present are never removed by a vertex added
//! elsewhere: a cavity of a new vertex does not cross them, only subfaces which the vertex replaces.
//! The tetrahedralization is therefore Delaunay everywhere but close to the boundary.
class RBoundaryRecovery
{

    public:

        //! Marker of a subface or a subsegment which does not belong to any group.
        static constexpr int noMarker = -1;

        enum VertexType
        {
            //! Vertex given on input.
            VertexInput = 0,
            //! Vertex added on a subsegment.
            VertexSegment,
            //! Vertex added on a subface.
            VertexFacet,
            //! Vertex added inside of a volume.
            VertexVolume,
            nVertexTypes
        };

        struct SubSegment
        {
            int v[2];
            int marker;
            bool alive;
            bool queued;
        };

        struct SubFace
        {
            int v[3];
            int marker;
            bool alive;
            bool queued;
        };

        enum EncroachmentType
        {
            EncroachmentNone = 0,
            EncroachmentSegment,
            EncroachmentFace
        };

        //! Subsegment or subface which a point lies too close to.
        struct Encroachment
        {
            EncroachmentType type = EncroachmentNone;
            int id = -1;
        };

    protected:

        struct FaceKey
        {
            int v[3];

            bool operator ==(const FaceKey &key) const
            {
                return this->v[0] == key.v[0] && this->v[1] == key.v[1] && this->v[2] == key.v[2];
            }
        };

        struct FaceKeyHash
        {
            size_t operator ()(const FaceKey &key) const
            {
                return size_t(  (uint64_t(key.v[0]) * 73856093ULL)
                              ^ (uint64_t(key.v[1]) * 19349663ULL)
                              ^ (uint64_t(key.v[2]) * 83492791ULL));
            }
        };

        enum FlipMode
        {
            //! Flip edge if its subfaces are not Delaunay.
            FlipDelaunay,
            //! Flip edge if vertices of its subfaces lie on one circle.
            FlipCocircular
        };

        //! Edge around subfaces which are replaced by subfaces connecting a new vertex.
        struct LoopEdge
        {
            int a;
            int b;
            //! Marker of the replaced subface.
            int marker;
        };

        //! Result of a search for a subface containing a point.
        struct FacetLocation
        {
            //! Subface which contains the point, -1 if a subsegment was hit before the point was reached.
            int faceID = -1;
            //! Subsegment which was hit or which the point lies on, -1 if there is none.
            int segmentID = -1;
            //! Whether the point lies on an edge of the subface.
            bool onEdge = false;
            //! Vertices of the edge which the point lies on.
            int edge[2] = {-1,-1};
            //! All visited subfaces.
            std::vector<int> visitedFaceIDs;
        };

        RTetraMesh &mesh;

        //! Maximum sine of an angle between normals of two subfaces lying in one plane.
        double tolerance;
        //! Length below which subsegments and subfaces are not split.
        double minLength;
        //! Maximum number of added vertices.
        uint64_t maxSteinerVertices;

        std::vector<SubSegment> segments;
        std::vector<SubFace> faces;

        //! Subsegment of each edge.
        std::unordered_map<uint64_t,int> segmentMap;
        //! Subface of each triple of vertices.
        std::unordered_map<FaceKey,int,FaceKeyHash> faceMap;
        //! Subfaces sharing each edge.
        std::unordered_map<uint64_t,std::vector<int>> edgeFaces;

        //! Subsegments which have to be checked for presence in the tetrahedralization.
        std::vector<int> segmentQueue;
        //! Subfaces which have to be checked for presence in the tetrahedralization.
        std::vector<int> faceQueue;

        //! Type of each vertex.
        std::vector<unsigned char> vertexTypes;
        //! Whether each vertex belongs to a subsegment or to a subface.
        std::vector<unsigned char> boundaryVertices;
        //! Required mesh size at each vertex, used only if positive.
        std::vector<double> vertexSizes;

        //! Marks of subfaces which are being replaced, such subfaces do not block a cavity.
        std::vector<unsigned int> replacedMarks;
        unsigned int replacedStamp;
        //! Marks of subfaces visited while a facet is searched.
        std::vector<unsigned int> facetMarks;
        unsigned int facetStamp;
        //! Function telling whether a face of a tetrahedron blocks a cavity.
        std::function<bool(int,int)> blockedFunction;

        //! Number of added vertices of each type.
        uint64_t nSteinerVertices[nVertexTypes];

        //! Whether created tetrahedra are collected.
        bool collectCreatedTetras;
        std::vector<int> createdTetras;

    private:

        //! Add subsegment, or set its marker if it exists. Return its position.
        int addSubSegment(int a, int b, int marker);

        //! Remove subsegment.
        void removeSubSegment(int segmentID);

        //! Add subface. Return its position, -1 if it is degenerated or if it exists.
        int addSubFace(int a, int b, int c, int marker);

        //! Remove subface.
        void removeSubFace(int faceID);

        //! Put subsegment to a queue of subsegments to be checked.
        void queueSubSegment(int segmentID);

        //! Put subface to a queue of subfaces to be checked.
        void queueSubFace(int faceID);

        //! Put subsegments and subfaces of the current cavity to queues.
        void queueCavity();

        //! Find mesh size at a point from vertices of a tetrahedron which contains it.
        double interpolateSize(const double *p, int tetraID) const;

        //! Return true if a face of a tetrahedron is a subface which is not being replaced.
        bool isFaceBlocked(int tetraID, int face) const;

        //! Find subfaces of a facet which have to be replaced to add a vertex: those which circumscribed
        //! circles contain it and which can be reached from given subfaces the vertex lies in.
        //! Edge with given key is not crossed. No subfaces are returned if the vertex can not be added.
        void findFacetCavity(const double *p, const std::vector<int> &startFaceIDs, uint64_t exceptKey, std::vector<int> &removedFaceIDs);

        //! Find edges around subfaces which are replaced, without an edge with given key.
        std::vector<LoopEdge> findLoopEdges(const std::vector<int> &removedFaceIDs, uint64_t exceptKey);

        //! Find cavity of the mesh for a vertex which replaces given subfaces and splits edge (a,b), if given.
        bool findBoundaryCavity(const double *p, const std::vector<int> &removedFaceIDs, int a, int b, int locatedTetraID);

        //! Return subsegment which would be removed if the current cavity of the mesh would be filled, -1 if there is none.
        int findLostSubSegment(int exceptSegmentID) const;

        //! Return distance of a point from the closest vertex of the current cavity of the mesh.
        double findCavityRadius(const double *p) const;

        //! Flip an edge shared by two subfaces if it is allowed by given mode.
        //! Vertices opposite to the edge are returned.
        bool flipEdge(int a, int b, FlipMode mode, int &c, int &d);

        //! Flip edges of a facet until its subfaces are Delaunay.
        void restoreDelaunay(std::vector<uint64_t> &edgeStack);

        //! Change triangulation of a facet to the one found in the tetrahedralization.
        bool matchSubFace(int faceID);

        //! Find subface of the same facet which contains a point lying in plane of given subface.
        FacetLocation locateInFacet(int faceID, const double *p) const;

        //! Replace each subface at an edge by two subfaces with a vertex on the edge.
        void splitEdgeFaces(int a, int b, int vertexID, std::vector<uint64_t> &edgeStack);

        //! Recover subface which is missing in the tetrahedralization.
        void recoverSubFace(int faceID);

    public:

        //! Constructor.
        //! All points of the mesh are taken as input vertices.
        RBoundaryRecovery(RTetraMesh &mesh, double tolerance);

        //! Set length below which nothing is split.
        void setMinLength(double minLength);

        //! Set maximum number of added vertices.
        void setMaxSteinerVertices(uint64_t maxSteinerVertices);

        //! Return number of added vertices of given type.
        uint64_t getNSteinerVertices(VertexType type) const;

        //! Return number of all added vertices.
        uint64_t getNSteinerVertices() const;

        //! Return required mesh size at a vertex.
        double getVertexSize(int vertexID) const;

        //! Set required mesh size at a vertex.
        void setVertexSize(int vertexID, double size);

        //! Return type of a vertex.
        VertexType getVertexType(int vertexID) const;

        //! Add triangle of the boundary.
        void addTriangle(int a, int b, int c, int marker);

        //! Add edge which has to be present in the mesh.
        void addEdge(int a, int b, int marker);

        //! Find subsegments and make facets Delaunay.
        //! Must be called after all triangles and edges were added.
        void initialize();

        //! Recover all missing subsegments and subfaces.
        void recover();

        //! Return subsegments, including removed ones.
        const std::vector<SubSegment> &getSubSegments() const;

        //! Return subfaces, including removed ones.
        const std::vector<SubFace> &getSubFaces() const;

        //! Return position of a subsegment, -1 if there is none.
        int findSubSegment(int a, int b) const;

        //! Return position of a subface, -1 if there is none.
        int findSubFace(int a, int b, int c) const;

        //! Return position of a subface which is a face of a tetrahedron, -1 if there is none.
        int findSubFace(const RTetraMesh::Face &face) const;

        //! Return length of a subsegment.
        double findSubSegmentLength(int segmentID) const;

        //! Return radius of a circumscribed circle of a subface.
        double findSubFaceRadius(int faceID) const;

        //! Find cavity of the mesh for a point which does not lie on the boundary. The cavity does not cross subfaces.
        //! Return false if there is none.
        bool prepareVolumeCavity(const double *p, int seedTetraID);

        //! Find subsegment or subface which would be removed or which would be too close to a point,
        //! if the point would be inserted to the current cavity of the mesh.
        Encroachment findEncroachment(const double *p) const;

        //! Insert point to the current cavity of the mesh and return its position.
        int insertCavityPoint(const double *p, VertexType type);

        //! Split subsegment by a new vertex. Return false if it could not be split,
        //! or if the new vertex would be closer than minRadius to another one.
        bool splitSubSegment(int segmentID, double minRadius = 0.0);

        //! Split subface by a new vertex in center of its circumscribed circle, or split a subsegment which
        //! lies in the way to the center. Return false if nothing could be split,
        //! or if the new vertex would be closer than minRadius to another one.
        bool splitSubFace(int faceID, double minRadius = 0.0);

        //! Set regions of tetrahedra: tetrahedra which can be reached from outside without crossing a subface
        //! have none, each other group of tetrahedra enclosed by subfaces is a region.
        //! Return number of regions.
        int classifyRegions();

        //! Set whether created tetrahedra are collected.
        void setCollectCreatedTetras(bool collectCreatedTetras);

        //! Return collected tetrahedra. Some of them may be removed or replaced already.
        std::vector<int> &getCreatedTetras();

};

#endif // REL_BOUNDARY_RECOVERY_H
