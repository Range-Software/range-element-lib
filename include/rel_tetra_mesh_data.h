#ifndef REL_TETRA_MESH_DATA_H
#define REL_TETRA_MESH_DATA_H

#include <vector>

//! Input and output of the tetrahedral mesh generator.
//!
//! As an input it describes a piecewise linear complex: nodes, triangles of the boundary and of
//! internal interfaces, edges which have to be present in the mesh and optionally an existing
//! tetrahedral mesh to be reconstructed.
//! As an output it holds the generated mesh. Input nodes keep their positions in the list of
//! nodes, nodes added by the generator (Steiner nodes) follow them.
class RTetraMeshData
{

    public:

        //! Marker of an item which does not belong to any group.
        static constexpr int noMarker = -1;

        //! Position of a tetrahedron which does not exist.
        static constexpr int noTetrahedron = -1;

        struct Node
        {
            //! Coordinates.
            double x = 0.0;
            double y = 0.0;
            double z = 0.0;
            //! Marker of a group of points the node belongs to.
            int marker = noMarker;
            //! Required size of mesh at the node, used only if it is positive.
            double size = 0.0;
        };

        struct Edge
        {
            int nodeIDs[2] = {0,0};
            //! Marker of a group of lines the edge belongs to.
            int marker = noMarker;
        };

        struct Triangle
        {
            int nodeIDs[3] = {0,0,0};
            //! Marker of a group of surfaces the triangle belongs to.
            int marker = noMarker;
        };

        struct Tetrahedron
        {
            //! Nodes are ordered so, that the fourth one lies on the side of the first three,
            //! which their normal (n2-n1)x(n3-n1) points to.
            int nodeIDs[4] = {0,0,0,0};
            //! Marker of a group of volumes (region) the tetrahedron belongs to.
            int marker = noMarker;
            //! Tetrahedra sharing a face opposite to each node.
            int neighborIDs[4] = {noTetrahedron,noTetrahedron,noTetrahedron,noTetrahedron};
        };

    protected:

        std::vector<Node> nodes;
        std::vector<Edge> edges;
        std::vector<Triangle> triangles;
        std::vector<Tetrahedron> tetrahedra;

    public:

        //! Constructor.
        RTetraMeshData();

        //! Remove all nodes, edges, triangles and tetrahedra.
        void clear();

        //! Return const reference to nodes.
        const std::vector<Node> &getNodes() const;

        //! Return reference to nodes.
        std::vector<Node> &getNodes();

        //! Return const reference to edges.
        const std::vector<Edge> &getEdges() const;

        //! Return reference to edges.
        std::vector<Edge> &getEdges();

        //! Return const reference to triangles.
        const std::vector<Triangle> &getTriangles() const;

        //! Return reference to triangles.
        std::vector<Triangle> &getTriangles();

        //! Return const reference to tetrahedra.
        const std::vector<Tetrahedron> &getTetrahedra() const;

        //! Return reference to tetrahedra.
        std::vector<Tetrahedron> &getTetrahedra();

        //! Add node and return its position.
        int addNode(double x, double y, double z, int marker = noMarker, double size = 0.0);

        //! Add edge and return its position.
        int addEdge(int node1, int node2, int marker = noMarker);

        //! Add triangle and return its position.
        int addTriangle(int node1, int node2, int node3, int marker = noMarker);

        //! Add tetrahedron and return its position.
        int addTetrahedron(int node1, int node2, int node3, int node4, int marker = noMarker);

        //! Return volume of a tetrahedron, which is positive if its nodes are ordered correctly.
        double findTetrahedronVolume(int tetrahedronID) const;

        //! Return area of a triangle.
        double findTriangleArea(int triangleID) const;

};

#endif // REL_TETRA_MESH_DATA_H
