#include <cmath>

#include "rel_tetra_mesh_data.h"

RTetraMeshData::RTetraMeshData()
{
}

void RTetraMeshData::clear()
{
    this->nodes.clear();
    this->edges.clear();
    this->triangles.clear();
    this->tetrahedra.clear();
}

const std::vector<RTetraMeshData::Node> &RTetraMeshData::getNodes() const
{
    return this->nodes;
}

std::vector<RTetraMeshData::Node> &RTetraMeshData::getNodes()
{
    return this->nodes;
}

const std::vector<RTetraMeshData::Edge> &RTetraMeshData::getEdges() const
{
    return this->edges;
}

std::vector<RTetraMeshData::Edge> &RTetraMeshData::getEdges()
{
    return this->edges;
}

const std::vector<RTetraMeshData::Triangle> &RTetraMeshData::getTriangles() const
{
    return this->triangles;
}

std::vector<RTetraMeshData::Triangle> &RTetraMeshData::getTriangles()
{
    return this->triangles;
}

const std::vector<RTetraMeshData::Tetrahedron> &RTetraMeshData::getTetrahedra() const
{
    return this->tetrahedra;
}

std::vector<RTetraMeshData::Tetrahedron> &RTetraMeshData::getTetrahedra()
{
    return this->tetrahedra;
}

int RTetraMeshData::addNode(double x, double y, double z, int marker, double size)
{
    Node node;
    node.x = x;
    node.y = y;
    node.z = z;
    node.marker = marker;
    node.size = size;
    this->nodes.push_back(node);
    return int(this->nodes.size()) - 1;
}

int RTetraMeshData::addEdge(int node1, int node2, int marker)
{
    Edge edge;
    edge.nodeIDs[0] = node1;
    edge.nodeIDs[1] = node2;
    edge.marker = marker;
    this->edges.push_back(edge);
    return int(this->edges.size()) - 1;
}

int RTetraMeshData::addTriangle(int node1, int node2, int node3, int marker)
{
    Triangle triangle;
    triangle.nodeIDs[0] = node1;
    triangle.nodeIDs[1] = node2;
    triangle.nodeIDs[2] = node3;
    triangle.marker = marker;
    this->triangles.push_back(triangle);
    return int(this->triangles.size()) - 1;
}

int RTetraMeshData::addTetrahedron(int node1, int node2, int node3, int node4, int marker)
{
    Tetrahedron tetrahedron;
    tetrahedron.nodeIDs[0] = node1;
    tetrahedron.nodeIDs[1] = node2;
    tetrahedron.nodeIDs[2] = node3;
    tetrahedron.nodeIDs[3] = node4;
    tetrahedron.marker = marker;
    this->tetrahedra.push_back(tetrahedron);
    return int(this->tetrahedra.size()) - 1;
}

double RTetraMeshData::findTetrahedronVolume(int tetrahedronID) const
{
    const Tetrahedron &rTetrahedron = this->tetrahedra[size_t(tetrahedronID)];
    const Node &a = this->nodes[size_t(rTetrahedron.nodeIDs[0])];
    const Node &b = this->nodes[size_t(rTetrahedron.nodeIDs[1])];
    const Node &c = this->nodes[size_t(rTetrahedron.nodeIDs[2])];
    const Node &d = this->nodes[size_t(rTetrahedron.nodeIDs[3])];

    const double abx = b.x - a.x, aby = b.y - a.y, abz = b.z - a.z;
    const double acx = c.x - a.x, acy = c.y - a.y, acz = c.z - a.z;
    const double adx = d.x - a.x, ady = d.y - a.y, adz = d.z - a.z;

    return (  (aby * acz - abz * acy) * adx
            + (abz * acx - abx * acz) * ady
            + (abx * acy - aby * acx) * adz) / 6.0;
}

double RTetraMeshData::findTriangleArea(int triangleID) const
{
    const Triangle &rTriangle = this->triangles[size_t(triangleID)];
    const Node &a = this->nodes[size_t(rTriangle.nodeIDs[0])];
    const Node &b = this->nodes[size_t(rTriangle.nodeIDs[1])];
    const Node &c = this->nodes[size_t(rTriangle.nodeIDs[2])];

    const double abx = b.x - a.x, aby = b.y - a.y, abz = b.z - a.z;
    const double acx = c.x - a.x, acy = c.y - a.y, acz = c.z - a.z;

    return 0.5 * std::hypot(aby * acz - abz * acy,
                            abz * acx - abx * acz,
                            abx * acy - aby * acx);
}
