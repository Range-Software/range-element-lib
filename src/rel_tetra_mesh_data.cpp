#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_map>

#include "rel_tetra_mesh_data.h"

const int RTetraMeshData::triangleEdgeNodes[3][2] = {{0, 1}, {1, 2}, {2, 0}};

const int RTetraMeshData::tetrahedronEdgeNodes[6][2] = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};

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

int RTetraMeshData::getElementOrder() const
{
    for (const Tetrahedron &rTetrahedron : this->tetrahedra)
    {
        if (rTetrahedron.midNodeIDs[0] != noNode)
        {
            return 2;
        }
    }
    for (const Triangle &rTriangle : this->triangles)
    {
        if (rTriangle.midNodeIDs[0] != noNode)
        {
            return 2;
        }
    }
    for (const Edge &rEdge : this->edges)
    {
        if (rEdge.midNodeID != noNode)
        {
            return 2;
        }
    }
    return 1;
}

int RTetraMeshData::createMidEdgeNodes()
{
    this->removeMidEdgeNodes();

    const size_t nCornerNodes = this->nodes.size();

    // Node of each edge. Nodes are numbered in the order their edges are met,
    // which does not depend on the order of items in the table.
    std::unordered_map<uint64_t, int> edgeNodes;
    edgeNodes.reserve(2 * this->tetrahedra.size() + 2 * this->triangles.size() + this->edges.size());

    auto findMidNode = [&](int node1, int node2) -> int
    {
        const uint64_t key = (uint64_t(uint32_t(std::min(node1, node2))) << 32) | uint64_t(uint32_t(std::max(node1, node2)));
        auto iter = edgeNodes.find(key);
        if (iter != edgeNodes.end())
        {
            return iter->second;
        }
        // References are not kept, the list of nodes grows.
        const Node a = this->nodes[size_t(node1)];
        const Node b = this->nodes[size_t(node2)];
        const int nodeID = this->addNode(a.x + 0.5 * (b.x - a.x),
                                         a.y + 0.5 * (b.y - a.y),
                                         a.z + 0.5 * (b.z - a.z),
                                         noMarker,
                                         a.size > 0.0 && b.size > 0.0 ? 0.5 * (a.size + b.size) : 0.0);
        edgeNodes.emplace(key, nodeID);
        return nodeID;
    };

    for (Tetrahedron &rTetrahedron : this->tetrahedra)
    {
        for (int i = 0; i < 6; i++)
        {
            rTetrahedron.midNodeIDs[i] = findMidNode(rTetrahedron.nodeIDs[tetrahedronEdgeNodes[i][0]],
                                                     rTetrahedron.nodeIDs[tetrahedronEdgeNodes[i][1]]);
        }
    }
    for (Triangle &rTriangle : this->triangles)
    {
        for (int i = 0; i < 3; i++)
        {
            rTriangle.midNodeIDs[i] = findMidNode(rTriangle.nodeIDs[triangleEdgeNodes[i][0]],
                                                  rTriangle.nodeIDs[triangleEdgeNodes[i][1]]);
        }
    }
    for (Edge &rEdge : this->edges)
    {
        rEdge.midNodeID = findMidNode(rEdge.nodeIDs[0], rEdge.nodeIDs[1]);
    }

    return int(this->nodes.size() - nCornerNodes);
}

int RTetraMeshData::removeMidEdgeNodes()
{
    // Nodes which are a mid-edge node of an element and a corner node of none.
    std::vector<unsigned char> remove(this->nodes.size(), 0);
    bool found = false;

    auto mark = [&](int &rNodeID)
    {
        if (rNodeID != noNode)
        {
            remove[size_t(rNodeID)] = 1;
            rNodeID = noNode;
            found = true;
        }
    };

    for (Tetrahedron &rTetrahedron : this->tetrahedra)
    {
        for (int &rNodeID : rTetrahedron.midNodeIDs)
        {
            mark(rNodeID);
        }
    }
    for (Triangle &rTriangle : this->triangles)
    {
        for (int &rNodeID : rTriangle.midNodeIDs)
        {
            mark(rNodeID);
        }
    }
    for (Edge &rEdge : this->edges)
    {
        mark(rEdge.midNodeID);
    }
    if (!found)
    {
        return 0;
    }

    for (const Tetrahedron &rTetrahedron : this->tetrahedra)
    {
        for (int nodeID : rTetrahedron.nodeIDs)
        {
            remove[size_t(nodeID)] = 0;
        }
    }
    for (const Triangle &rTriangle : this->triangles)
    {
        for (int nodeID : rTriangle.nodeIDs)
        {
            remove[size_t(nodeID)] = 0;
        }
    }
    for (const Edge &rEdge : this->edges)
    {
        for (int nodeID : rEdge.nodeIDs)
        {
            remove[size_t(nodeID)] = 0;
        }
    }

    // Nodes which stay keep their order.
    std::vector<int> nodeMap(this->nodes.size(), noNode);
    size_t nKept = 0;
    for (size_t i = 0; i < this->nodes.size(); i++)
    {
        if (!remove[i])
        {
            nodeMap[i] = int(nKept);
            this->nodes[nKept++] = this->nodes[i];
        }
    }
    const int nRemoved = int(this->nodes.size() - nKept);
    this->nodes.resize(nKept);

    for (Tetrahedron &rTetrahedron : this->tetrahedra)
    {
        for (int &rNodeID : rTetrahedron.nodeIDs)
        {
            rNodeID = nodeMap[size_t(rNodeID)];
        }
    }
    for (Triangle &rTriangle : this->triangles)
    {
        for (int &rNodeID : rTriangle.nodeIDs)
        {
            rNodeID = nodeMap[size_t(rNodeID)];
        }
    }
    for (Edge &rEdge : this->edges)
    {
        for (int &rNodeID : rEdge.nodeIDs)
        {
            rNodeID = nodeMap[size_t(rNodeID)];
        }
    }
    return nRemoved;
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

    return ((aby * acz - abz * acy) * adx + (abz * acx - abx * acz) * ady + (abx * acy - aby * acx) * adz) / 6.0;
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
