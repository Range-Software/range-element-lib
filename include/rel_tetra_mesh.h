#ifndef REL_TETRA_MESH_H
#define REL_TETRA_MESH_H

#include <cstdint>
#include <functional>
#include <vector>

//! Tetrahedralization of a set of points which points can be added to.
//!
//! It is a Delaunay tetrahedralization as long as points are inserted without restrictions
//! (see findCavity()), ties of points lying on one sphere are broken by symbolic perturbation.
//!
//! Space outside of the convex hull of points is covered by ghost tetrahedra, each made of a
//! face of the hull and of a ghost vertex lying in infinity. Each tetrahedron has therefore
//! four neighbors and no special treatment of the hull is needed.
class RTetraMesh
{

    public:

        //! Vertex lying in infinity.
        static constexpr int ghostVertex = -1;
        //! First vertex of a tetrahedron which was removed.
        static constexpr int deadVertex = -2;
        //! Tetrahedron which does not exist.
        static constexpr int noTetra = -1;
        //! Region of a tetrahedron lying outside of meshed domain.
        static constexpr int noRegion = -1;

        struct Point
        {
            double x[3];
        };

        struct Tetra
        {
            //! Vertices. Vertex v[3] lies on the side of (v[0],v[1],v[2]) which its normal points to.
            int v[4];
            //! Neighbors, n[i] shares the face opposite to v[i].
            int n[4];
            //! Region which the tetrahedron belongs to.
            int region;
        };

        //! Face of a tetrahedron given by a position of its opposite vertex.
        struct Face
        {
            int tetra;
            int face;
        };

        //! Positions of vertices of a face opposite to each vertex of a tetrahedron.
        //! They are ordered so that the opposite vertex lies on the side which their normal points to.
        static const int faceVertices[4][3];

        //! Positions of vertices of edges of a tetrahedron.
        static const int edgeVertices[6][2];

    protected:

        std::vector<Point> points;
        std::vector<Tetra> tetras;
        //! Positions of removed tetrahedra which can be used again.
        std::vector<int> freeTetras;
        //! Tetrahedron each vertex belongs to.
        std::vector<int> vertexTetras;
        //! Number of tetrahedra which were not removed, including ghost ones.
        uint64_t nTetras;
        //! Tetrahedron which a search for a point starts from by default.
        int lastTetra;

        //! Marks of tetrahedra visited while a cavity is searched.
        std::vector<unsigned int> cavityMarks;
        unsigned int cavityStamp;
        //! Marks of tetrahedra visited while a star of a vertex is searched.
        mutable std::vector<unsigned int> starMarks;
        mutable unsigned int starStamp;
        mutable std::vector<int> starBuffer;
        //! State of a random number generator used by a search for a point.
        mutable unsigned int walkSeed;
        //! Relative distance of a point from a plane, below which the point is considered to lie in the plane.
        double coplanarTolerance;

        //! Tetrahedra which have to be removed to add a point.
        std::vector<int> cavityTetras;
        //! Faces of cavity tetrahedra which their neighbors do not belong to the cavity.
        std::vector<Face> cavityFaces;
        //! Edge of a cavity boundary face, which a new tetrahedron is connected to another one through.
        struct FaceLink
        {
            //! Key of the edge.
            uint64_t key;
            //! Position of the boundary face.
            int face;
            //! Position of the neighbor in the new tetrahedron.
            int side;
        };
        //! Edges of boundary faces of the cavity sorted by their keys.
        std::vector<FaceLink> cavityLinks;
        //! Tetrahedra created when the last cavity was filled.
        std::vector<int> createdTetras;

    private:

        //! Return position of a new tetrahedron.
        int allocateTetra();

    public:

        //! Constructor.
        RTetraMesh();

        //! Remove all points and tetrahedra.
        void clear();

        //! Add point without inserting it into the tetrahedralization and return its position.
        int addPoint(double x, double y, double z);

        //! Return number of points.
        int getNPoints() const;

        //! Return coordinates of a point.
        inline const double *getPoint(int pointID) const
        {
            return this->points[size_t(pointID)].x;
        }

        //! Return number of positions of tetrahedra, including removed ones.
        int getNTetraSlots() const;

        //! Return number of tetrahedra, including ghost ones.
        uint64_t getNTetras() const;

        //! Return tetrahedron.
        inline const Tetra &getTetra(int tetraID) const
        {
            return this->tetras[size_t(tetraID)];
        }

        //! Return true if tetrahedron was not removed.
        inline bool isAlive(int tetraID) const
        {
            return this->tetras[size_t(tetraID)].v[0] != deadVertex;
        }

        //! Return position of the ghost vertex in a tetrahedron or -1 if there is none.
        inline int findGhostIndex(int tetraID) const
        {
            const Tetra &rTetra = this->tetras[size_t(tetraID)];
            for (int i=0;i<4;i++)
            {
                if (rTetra.v[i] == ghostVertex)
                {
                    return i;
                }
            }
            return -1;
        }

        //! Return true if tetrahedron has the ghost vertex.
        inline bool isGhost(int tetraID) const
        {
            return this->findGhostIndex(tetraID) >= 0;
        }

        //! Set region of a tetrahedron.
        void setRegion(int tetraID, int region);

        //! Return tetrahedron a vertex belongs to, noTetra if it was not inserted yet.
        int getVertexTetra(int pointID) const;

        //! Create the first tetrahedron and ghost tetrahedra around it.
        //! Points must not be coplanar.
        void createFirstTetra(int a, int b, int c, int d);

        //! Find tetrahedron which contains a point (inside or on its boundary), starting at given one.
        //! If point lies outside of the convex hull, ghost tetrahedron which it conflicts with is returned.
        int locate(const double *p, int startTetraID = noTetra) const;

        //! Find vertex of a tetrahedron with the same coordinates as a point, return -1 if there is none.
        int findDuplicateVertex(int tetraID, const double *p) const;

        //! Set relative distance of a point from a plane, below which the point is considered to lie in the plane.
        //! It is used where a decision based on an exact position would produce flat tetrahedra,
        //! e.g. from points of a planar face which coordinates are rounded.
        void setCoplanarTolerance(double coplanarTolerance);

        //! Return orientation of four points (see RGeometricPredicates::orient3d),
        //! zero if the fourth point lies in plane of the first three or closer to it than the tolerance.
        double findOrientation(const double *pa, const double *pb, const double *pc, const double *pd) const;

        //! Return true if point lies inside of the circumscribed sphere of a tetrahedron.
        bool isInConflict(int tetraID, const double *p) const;

        //! Find tetrahedra which have to be removed to add a point (cavity).
        //! The search starts at given tetrahedron and goes through tetrahedra which the point conflicts with.
        //! Forced tetrahedra are always a part of the cavity. Faces for which isBlocked(tetraID,face) returns
        //! true are not crossed and never end inside of the cavity.
        //! Tetrahedra are then removed from the cavity until the point can see its whole boundary,
        //! so that all new tetrahedra have a positive volume.
        //! Return false if there is no such cavity. Then the mesh can not be changed by fillCavity().
        bool findCavity(const double *p,
                        int seedTetraID,
                        const std::vector<int> *pForcedTetras = nullptr,
                        const std::function<bool(int,int)> *pIsBlocked = nullptr);

        //! Return tetrahedra of the last found cavity.
        const std::vector<int> &getCavityTetras() const;

        //! Return boundary faces of the last found cavity.
        const std::vector<Face> &getCavityFaces() const;

        //! Return true if tetrahedron belongs to the last found cavity.
        inline bool isInCavity(int tetraID) const
        {
            const unsigned int mark = this->cavityMarks[size_t(tetraID)];
            return mark == this->cavityStamp || mark == this->cavityStamp + 2;
        }

        //! Replace tetrahedra of the last found cavity by tetrahedra connecting its boundary to a point.
        //! New tetrahedra take region of the removed tetrahedron they share a face with.
        void fillCavity(int pointID);

        //! Return tetrahedra created when the last cavity was filled.
        const std::vector<int> &getCreatedTetras() const;

        //! Insert point into the tetrahedralization.
        //! Return false if a vertex with the same coordinates exists, or if the point is too close to
        //! a vertex to form tetrahedra with it. Position of that vertex is then set to duplicateID.
        bool insertVertex(int pointID, int hintTetraID, int &duplicateID);

        //! Find all tetrahedra sharing a vertex.
        void findStar(int pointID, std::vector<int> &star) const;

        //! Find tetrahedron having an edge, return noTetra if there is none.
        int findEdge(int a, int b) const;

        //! Find tetrahedron having a face, return false if there is none.
        bool findFace(int a, int b, int c, Face &face) const;

};

#endif // REL_TETRA_MESH_H
