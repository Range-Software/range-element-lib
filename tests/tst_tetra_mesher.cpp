#include <QTest>
#include <QFile>
#include <QTemporaryDir>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QRegularExpression>
#include <QTextStream>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <vector>

#include <rbl_error.h>

#include "rel_tetra_mesh_data.h"
#include "rel_tetra_mesh_options.h"
#include "rel_tetra_mesher.h"
#include "rel_tetra_mesh_quality.h"

class TestTetraMesher : public QObject
{
    Q_OBJECT

    private:

        //! Add surface of a box. Normals of triangles point out of the box.
        //! If flipSide is set, the side at x0 is divided by the other diagonal, the same as the side at x1 is.
        static void addBox(RTetraMeshData &data, double x0, double y0, double z0, double x1, double y1, double z1, int marker, bool flipSide = false)
        {
            const int n = int(data.getNodes().size());
            data.addNode(x0,y0,z0);
            data.addNode(x1,y0,z0);
            data.addNode(x1,y1,z0);
            data.addNode(x0,y1,z0);
            data.addNode(x0,y0,z1);
            data.addNode(x1,y0,z1);
            data.addNode(x1,y1,z1);
            data.addNode(x0,y1,z1);
            int triangles[12][3] = {{0,2,1},{0,3,2},{4,5,6},{4,6,7},
                                    {0,1,5},{0,5,4},{1,2,6},{1,6,5},
                                    {2,3,7},{2,7,6},{3,0,4},{3,4,7}};
            if (flipSide)
            {
                const int flipped[2][3] = {{0,4,7},{0,7,3}};
                for (int i=0;i<2;i++)
                {
                    for (int j=0;j<3;j++)
                    {
                        triangles[10+i][j] = flipped[i][j];
                    }
                }
            }
            for (const auto &triangle : triangles)
            {
                data.addTriangle(n+triangle[0],n+triangle[1],n+triangle[2],marker);
            }
        }

        //! Add surface of a box 2 x 1 x 1 made of two cubes, and a wall between the cubes if required.
        //! The surface has marker 0, the wall has marker 1.
        static void addDividedBox(RTetraMeshData &data, bool addWall)
        {
            auto node = [](int i, int j, int k) { return i*4 + j*2 + k; };
            for (int i=0;i<3;i++)
            {
                for (int j=0;j<2;j++)
                {
                    for (int k=0;k<2;k++)
                    {
                        data.addNode(i,j,k);
                    }
                }
            }
            auto addQuad = [&data](int a, int b, int c, int d, int marker)
            {
                data.addTriangle(a,b,c,marker);
                data.addTriangle(a,c,d,marker);
            };
            for (int i=0;i<2;i++)
            {
                for (int s=0;s<2;s++)
                {
                    addQuad(node(i,0,s),node(i+1,0,s),node(i+1,1,s),node(i,1,s),0);
                    addQuad(node(i,s,0),node(i+1,s,0),node(i+1,s,1),node(i,s,1),0);
                }
            }
            for (int i=0;i<3;i++)
            {
                if (i != 1 || addWall)
                {
                    addQuad(node(i,0,0),node(i,1,0),node(i,1,1),node(i,0,1),i == 1 ? 1 : 0);
                }
            }
        }

        //! Add surface of a prism which base is a polygon in plane XY.
        static void addPrism(RTetraMeshData &data, const std::vector<std::array<double,2>> &polygon, const std::vector<std::array<int,3>> &baseTriangles, double height)
        {
            const int n = int(polygon.size());
            for (int k=0;k<2;k++)
            {
                for (const auto &point : polygon)
                {
                    data.addNode(point[0],point[1],k*height);
                }
            }
            for (const auto &triangle : baseTriangles)
            {
                data.addTriangle(triangle[0],triangle[2],triangle[1],0);
                data.addTriangle(n+triangle[0],n+triangle[1],n+triangle[2],1);
            }
            for (int i=0;i<n;i++)
            {
                const int j = (i+1)%n;
                data.addTriangle(i,j,n+j,2);
                data.addTriangle(i,n+j,n+i,2);
            }
        }

        //! Add surface of a sphere with radius 1 made by dividing faces of an icosahedron.
        static void addSphere(RTetraMeshData &data, int nDivisions)
        {
            const double t = (1.0 + std::sqrt(5.0)) / 2.0;
            std::vector<std::array<double,3>> points = {{{-1,t,0}},{{1,t,0}},{{-1,-t,0}},{{1,-t,0}},{{0,-1,t}},{{0,1,t}},
                                                        {{0,-1,-t}},{{0,1,-t}},{{t,0,-1}},{{t,0,1}},{{-t,0,-1}},{{-t,0,1}}};
            std::vector<std::array<int,3>> faces = {{{0,11,5}},{{0,5,1}},{{0,1,7}},{{0,7,10}},{{0,10,11}},{{1,5,9}},{{5,11,4}},
                                                    {{11,10,2}},{{10,7,6}},{{7,1,8}},{{3,9,4}},{{3,4,2}},{{3,2,6}},{{3,6,8}},
                                                    {{3,8,9}},{{4,9,5}},{{2,4,11}},{{6,2,10}},{{8,6,7}},{{9,8,1}}};
            auto normalize = [](std::array<double,3> &point)
            {
                const double length = std::hypot(point[0],point[1],point[2]);
                for (double &rValue : point)
                {
                    rValue /= length;
                }
            };
            for (auto &point : points)
            {
                normalize(point);
            }
            for (int division=0;division<nDivisions;division++)
            {
                std::map<std::array<int,2>,int> middles;
                auto findMiddle = [&](int a, int b) -> int
                {
                    const std::array<int,2> key = {std::min(a,b),std::max(a,b)};
                    auto iter = middles.find(key);
                    if (iter != middles.end())
                    {
                        return iter->second;
                    }
                    std::array<double,3> point = {points[size_t(a)][0] + points[size_t(b)][0],
                                                  points[size_t(a)][1] + points[size_t(b)][1],
                                                  points[size_t(a)][2] + points[size_t(b)][2]};
                    normalize(point);
                    points.push_back(point);
                    middles[key] = int(points.size()) - 1;
                    return int(points.size()) - 1;
                };
                std::vector<std::array<int,3>> divided;
                for (const auto &face : faces)
                {
                    const int a = findMiddle(face[0],face[1]);
                    const int b = findMiddle(face[1],face[2]);
                    const int c = findMiddle(face[2],face[0]);
                    divided.push_back({{face[0],a,c}});
                    divided.push_back({{face[1],b,a}});
                    divided.push_back({{face[2],c,b}});
                    divided.push_back({{a,b,c}});
                }
                faces = divided;
            }
            const int n = int(data.getNodes().size());
            for (const auto &point : points)
            {
                data.addNode(point[0],point[1],point[2]);
            }
            for (const auto &face : faces)
            {
                data.addTriangle(n+face[0],n+face[1],n+face[2],0);
            }
        }

        //! Return the smallest angle between faces of a tetrahedron in degrees.
        static double findMinDihedralAngle(const RTetraMeshData &data, int tetrahedronID)
        {
            static const int faceNodes[4][3] = {{1,3,2},{0,2,3},{0,3,1},{0,1,2}};
            const RTetraMeshData::Tetrahedron &rTetrahedron = data.getTetrahedra()[size_t(tetrahedronID)];
            double normals[4][3];
            for (int i=0;i<4;i++)
            {
                const RTetraMeshData::Node &a = data.getNodes()[size_t(rTetrahedron.nodeIDs[faceNodes[i][0]])];
                const RTetraMeshData::Node &b = data.getNodes()[size_t(rTetrahedron.nodeIDs[faceNodes[i][1]])];
                const RTetraMeshData::Node &c = data.getNodes()[size_t(rTetrahedron.nodeIDs[faceNodes[i][2]])];
                const double u[3] = {b.x-a.x,b.y-a.y,b.z-a.z};
                const double v[3] = {c.x-a.x,c.y-a.y,c.z-a.z};
                normals[i][0] = u[1]*v[2] - u[2]*v[1];
                normals[i][1] = u[2]*v[0] - u[0]*v[2];
                normals[i][2] = u[0]*v[1] - u[1]*v[0];
                const double length = std::hypot(normals[i][0],normals[i][1],normals[i][2]);
                for (double &rValue : normals[i])
                {
                    rValue /= length;
                }
            }
            double minAngle = 180.0;
            for (int i=0;i<4;i++)
            {
                for (int j=i+1;j<4;j++)
                {
                    const double angleCos = -(normals[i][0]*normals[j][0] + normals[i][1]*normals[j][1] + normals[i][2]*normals[j][2]);
                    minAngle = std::min(minAngle,std::acos(std::max(-1.0,std::min(1.0,angleCos))) * 180.0 / 3.14159265358979323846);
                }
            }
            return minAngle;
        }

        static double findTotalVolume(const RTetraMeshData &data, int marker = -2)
        {
            double volume = 0.0;
            for (int i=0;i<int(data.getTetrahedra().size());i++)
            {
                if (marker == -2 || data.getTetrahedra()[size_t(i)].marker == marker)
                {
                    volume += data.findTetrahedronVolume(i);
                }
            }
            return volume;
        }

        static double findTotalArea(const RTetraMeshData &data, int marker = -2)
        {
            double area = 0.0;
            for (int i=0;i<int(data.getTriangles().size());i++)
            {
                if (marker == -2 || data.getTriangles()[size_t(i)].marker == marker)
                {
                    area += data.findTriangleArea(i);
                }
            }
            return area;
        }

        static double findDistance(const RTetraMeshData &data, int node1, int node2)
        {
            const RTetraMeshData::Node &a = data.getNodes()[size_t(node1)];
            const RTetraMeshData::Node &b = data.getNodes()[size_t(node2)];
            return std::hypot(a.x-b.x,a.y-b.y,a.z-b.z);
        }

        //! Return ratio of circumscribed sphere radius to the shortest edge of a tetrahedron.
        static double findRadiusEdgeRatio(const RTetraMeshData &data, int tetrahedronID)
        {
            const RTetraMeshData::Tetrahedron &rTetrahedron = data.getTetrahedra()[size_t(tetrahedronID)];
            const RTetraMeshData::Node &a = data.getNodes()[size_t(rTetrahedron.nodeIDs[0])];
            double e[3][3];
            for (int i=0;i<3;i++)
            {
                const RTetraMeshData::Node &b = data.getNodes()[size_t(rTetrahedron.nodeIDs[i+1])];
                e[i][0] = b.x - a.x;
                e[i][1] = b.y - a.y;
                e[i][2] = b.z - a.z;
            }
            auto cross = [](const double *u, const double *v, double *r)
            {
                r[0] = u[1]*v[2] - u[2]*v[1];
                r[1] = u[2]*v[0] - u[0]*v[2];
                r[2] = u[0]*v[1] - u[1]*v[0];
            };
            auto dot = [](const double *u, const double *v) { return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]; };
            double vw[3], wu[3], uv[3];
            cross(e[1],e[2],vw);
            cross(e[2],e[0],wu);
            cross(e[0],e[1],uv);
            const double det = dot(e[0],vw);
            double r[3];
            for (int i=0;i<3;i++)
            {
                r[i] = (dot(e[0],e[0])*vw[i] + dot(e[1],e[1])*wu[i] + dot(e[2],e[2])*uv[i]) / (2.0*det);
            }
            double minEdge = 1.0e300;
            for (int i=0;i<4;i++)
            {
                for (int j=i+1;j<4;j++)
                {
                    minEdge = std::min(minEdge,findDistance(data,rTetrahedron.nodeIDs[i],rTetrahedron.nodeIDs[j]));
                }
            }
            return std::sqrt(dot(r,r)) / minEdge;
        }

        //! Check that a mesh is valid: tetrahedra have positive volume, neighbors refer to each other,
        //! each face without a neighbor is one of triangles and input nodes were not changed.
        static void verifyMesh(const RTetraMeshData &input, const RTetraMeshData &output)
        {
            static const int faceNodes[4][3] = {{1,3,2},{0,2,3},{0,3,1},{0,1,2}};

            QVERIFY(output.getNodes().size() >= input.getNodes().size());
            for (size_t i=0;i<input.getNodes().size();i++)
            {
                QCOMPARE(output.getNodes()[i].x,input.getNodes()[i].x);
                QCOMPARE(output.getNodes()[i].y,input.getNodes()[i].y);
                QCOMPARE(output.getNodes()[i].z,input.getNodes()[i].z);
                QCOMPARE(output.getNodes()[i].marker,input.getNodes()[i].marker);
            }

            std::set<std::array<int,3>> triangles;
            for (const RTetraMeshData::Triangle &rTriangle : output.getTriangles())
            {
                std::array<int,3> key = {rTriangle.nodeIDs[0],rTriangle.nodeIDs[1],rTriangle.nodeIDs[2]};
                std::sort(key.begin(),key.end());
                // Each triangle is there only once.
                QVERIFY(triangles.insert(key).second);
            }

            const int nTetrahedra = int(output.getTetrahedra().size());
            QVERIFY(nTetrahedra > 0);

            std::set<std::array<int,3>> usedTriangles;
            for (int i=0;i<nTetrahedra;i++)
            {
                const RTetraMeshData::Tetrahedron &rTetrahedron = output.getTetrahedra()[size_t(i)];
                QVERIFY(output.findTetrahedronVolume(i) > 0.0);
                for (int j=0;j<4;j++)
                {
                    std::array<int,3> key = {rTetrahedron.nodeIDs[faceNodes[j][0]],
                                             rTetrahedron.nodeIDs[faceNodes[j][1]],
                                             rTetrahedron.nodeIDs[faceNodes[j][2]]};
                    std::sort(key.begin(),key.end());

                    const int neighborID = rTetrahedron.neighborIDs[j];
                    if (neighborID == RTetraMeshData::noTetrahedron)
                    {
                        QVERIFY(triangles.find(key) != triangles.end());
                        usedTriangles.insert(key);
                        continue;
                    }
                    QVERIFY(neighborID >= 0 && neighborID < nTetrahedra);
                    const RTetraMeshData::Tetrahedron &rNeighbor = output.getTetrahedra()[size_t(neighborID)];
                    bool found = false;
                    for (int k=0;k<4;k++)
                    {
                        if (rNeighbor.neighborIDs[k] != i)
                        {
                            continue;
                        }
                        std::array<int,3> neighborKey = {rNeighbor.nodeIDs[faceNodes[k][0]],
                                                         rNeighbor.nodeIDs[faceNodes[k][1]],
                                                         rNeighbor.nodeIDs[faceNodes[k][2]]};
                        std::sort(neighborKey.begin(),neighborKey.end());
                        found = found || neighborKey == key;
                    }
                    QVERIFY(found);
                    // Regions are separated by triangles.
                    if (rNeighbor.marker != rTetrahedron.marker)
                    {
                        QVERIFY(triangles.find(key) != triangles.end());
                    }
                    if (triangles.find(key) != triangles.end())
                    {
                        usedTriangles.insert(key);
                    }
                }
            }
            // Each triangle is a face of a tetrahedron.
            QCOMPARE(usedTriangles.size(),triangles.size());
        }

        //! Generate mesh and print an error if it fails.
        static bool generate(const RTetraMeshData &input, const RTetraMeshOptions &options, RTetraMeshData &output)
        {
            try
            {
                RTetraMesher::generate(input,options,output);
            }
            catch (const RError &rError)
            {
                qWarning("%s",rError.getMessage().toUtf8().constData());
                return false;
            }
            return true;
        }

        static RTetraMeshOptions createOptions(bool qualityMesh)
        {
            RTetraMeshOptions options;
            options.setVerbose(false);
            options.setQualityMesh(qualityMesh);
            return options;
        }

private slots:

    void chabot_data();
    void chabot();
    void box();
    void boxQuality();
    void boxSizeFunction();
    void concavePrism();
    void sharpWedge();
    void octahedron();
    void sphere();
    void sphereMinDihedralAngle();
    void twoRegions();
    void nestedBoxes();
    void pointsAndEdges();
    void quadrilateralFaces();
    void reconstruct();
    void duplicateNodes();
    void invalidInput();
};

void TestTetraMesher::chabot_data()
{
    QTest::addColumn<bool>("qualityMesh");
    QTest::addColumn<bool>("applicationSettings");
    QTest::newRow("boundary") << false << false;
    QTest::newRow("quality") << true << false;
    QTest::newRow("application") << true << true;
}

void TestTetraMesher::chabot()
{
    QFETCH(bool,qualityMesh);
    QFETCH(bool,applicationSettings);
    QFile file(QFINDTESTDATA("data/Chabot.stl"));
    QVERIFY(file.open(QIODevice::ReadOnly | QIODevice::Text));
    QTextStream stream(&file);
    RTetraMeshData input, output;
    std::map<std::array<double,3>,int> nodes;
    std::vector<int> triangle;
    while (!stream.atEnd())
    {
        const QString line = stream.readLine().trimmed();
        if (!line.startsWith("vertex "))
        {
            continue;
        }
        const QStringList fields = line.split(QRegularExpression("\\s+"),Qt::SkipEmptyParts);
        const std::array<double,3> point = {fields[1].toDouble(),fields[2].toDouble(),fields[3].toDouble()};
        auto entry = nodes.emplace(point,int(nodes.size()));
        if (entry.second)
        {
            input.addNode(point[0],point[1],point[2]);
        }
        triangle.push_back(entry.first->second);
        if (triangle.size() == 3)
        {
            input.addTriangle(triangle[0],triangle[1],triangle[2],0);
            triangle.clear();
        }
    }
    QCOMPARE(input.getNodes().size(),size_t(372));
    QCOMPARE(input.getTriangles().size(),size_t(756));
    auto options = createOptions(qualityMesh);
    if (applicationSettings)
    {
        options.setRadiusEdgeRatio(1.6);
        options.setMinDihedralAngle(10.0);
        options.setTolerance(1.0e-10);
        options.setVolumeConstraint(0.1);
    }
    QTemporaryDir diagnostics;
    if (applicationSettings)
    {
        QVERIFY(diagnostics.isValid());
        options.setQualityReportPath(diagnostics.filePath("quality.json").toStdString());
    }
    // Leave the default node budget: refinement must finish on its own.
    QVERIFY(generate(input,options,output));
    verifyMesh(input,output);
    qInfo("Chabot: %zu nodes, %zu tetrahedra",output.getNodes().size(),output.getTetrahedra().size());
    QVERIFY(output.getNodes().size() < size_t(10000));
    if (applicationSettings)
    {
        QFile report(QString::fromStdString(options.getQualityReportPath()));QVERIFY(report.open(QIODevice::ReadOnly));
        const auto document=QJsonDocument::fromJson(report.readAll()).object();
        const auto profile=document["profile"].toObject();
        QVERIFY(profile["total"].toDouble()>0);QVERIFY(profile["refinement"].toDouble()>0);
        QVERIFY(profile["qualityCacheHits"].toDouble()>0);QVERIFY(profile["starCacheHits"].toDouble()>0);
        QVERIFY(profile["trialReplacements"].toDouble()>0);QVERIFY(profile["journalBytes"].toDouble()>0);
        const auto patches=document["patches"].toArray();
        QVERIFY(!patches.empty());QVERIFY(patches.size()<=50);
        const auto quality=RTetraMeshQuality::summarize(output,options);
        QCOMPARE(quality.invalid,size_t(0));
        QCOMPARE(quality.nearFlat,size_t(0));
        QVERIFY(quality.maxRadiusEdge<1e4);
        QVERIFY(quality.radiusEdge95<3.5);
        QVERIFY(quality.badAngle*100<quality.elements*2);
        QVERIFY(quality.freeBelowOne*1000<quality.elements*3);
        QVERIFY(quality.minAngle>1e-3);
    }
    QVERIFY(std::fabs(findTotalArea(output) - findTotalArea(input)) < 1.0e-10);
    double volume = 0.0;
    for (const auto &face : input.getTriangles())
    {
        const auto &a = input.getNodes()[size_t(face.nodeIDs[0])];
        const auto &b = input.getNodes()[size_t(face.nodeIDs[1])];
        const auto &c = input.getNodes()[size_t(face.nodeIDs[2])];
        volume += (a.x*(b.y*c.z-b.z*c.y) + a.y*(b.z*c.x-b.x*c.z) + a.z*(b.x*c.y-b.y*c.x)) / 6.0;
    }
    QVERIFY(std::fabs(findTotalVolume(output) - std::fabs(volume)) < 1.0e-12);
}

void TestTetraMesher::box()
{
    RTetraMeshData input, output;
    addBox(input,0.0,0.0,0.0,1.0,1.0,1.0,3);

    QVERIFY(generate(input,createOptions(false),output));

    verifyMesh(input,output);
    QVERIFY(std::fabs(findTotalVolume(output) - 1.0) < 1.0e-12);
    QVERIFY(std::fabs(findTotalArea(output) - 6.0) < 1.0e-12);
    QVERIFY(std::fabs(findTotalArea(output,3) - 6.0) < 1.0e-12);
    for (const RTetraMeshData::Tetrahedron &rTetrahedron : output.getTetrahedra())
    {
        QCOMPARE(rTetrahedron.marker,0);
    }
    // Faces of a convex shape are a part of Delaunay tetrahedralization, no nodes are needed.
    QCOMPARE(output.getNodes().size(),input.getNodes().size());
    QVERIFY(output.getEdges().empty());
}

void TestTetraMesher::boxQuality()
{
    RTetraMeshData input, output;
    addBox(input,0.0,0.0,0.0,2.0,1.0,1.0,0);

    const double maxVolume = 0.005;
    const double maxRatio = 2.0;

    RTetraMeshOptions options = createOptions(true);
    options.setVolumeConstraint(maxVolume);
    options.setRadiusEdgeRatio(maxRatio);

    QVERIFY(generate(input,options,output));

    verifyMesh(input,output);
    QVERIFY(std::fabs(findTotalVolume(output) - 2.0) < 1.0e-10);
    QVERIFY(std::fabs(findTotalArea(output) - 10.0) < 1.0e-10);
    QVERIFY(output.getTetrahedra().size() >= size_t(2.0 / maxVolume));

    for (int i=0;i<int(output.getTetrahedra().size());i++)
    {
        QVERIFY(output.findTetrahedronVolume(i) <= maxVolume * (1.0 + 1.0e-12));
        QVERIFY(findRadiusEdgeRatio(output,i) <= maxRatio * (1.0 + 1.0e-12));
    }
}

void TestTetraMesher::boxSizeFunction()
{
    RTetraMeshData input, coarse, fine;
    addBox(input,0.0,0.0,0.0,1.0,1.0,1.0,0);

    RTetraMeshOptions options = createOptions(true);
    QVERIFY(generate(input,options,coarse));

    const double size = 0.15;
    for (RTetraMeshData::Node &rNode : input.getNodes())
    {
        rNode.size = size;
    }
    options.setUseSizeFunction(true);
    QVERIFY(generate(input,options,fine));

    verifyMesh(input,fine);
    QVERIFY(std::fabs(findTotalVolume(fine) - 1.0) < 1.0e-10);
    QVERIFY(fine.getTetrahedra().size() > 10 * coarse.getTetrahedra().size());

    for (const RTetraMeshData::Tetrahedron &rTetrahedron : fine.getTetrahedra())
    {
        for (int i=0;i<4;i++)
        {
            for (int j=i+1;j<4;j++)
            {
                QVERIFY(findDistance(fine,rTetrahedron.nodeIDs[i],rTetrahedron.nodeIDs[j]) <= 2.0 * size);
            }
        }
    }
    // Size is passed to added nodes.
    for (const RTetraMeshData::Node &rNode : fine.getNodes())
    {
        QVERIFY(std::fabs(rNode.size - size) < 1.0e-12);
    }
}

void TestTetraMesher::concavePrism()
{
    // Prism with a base of the shape of letter L, tetrahedra must not fill its inner corner.
    RTetraMeshData input, output;
    addPrism(input,
             {{{0.0,0.0}},{{2.0,0.0}},{{2.0,1.0}},{{1.0,1.0}},{{1.0,2.0}},{{0.0,2.0}}},
             {{{0,1,2}},{{0,2,3}},{{0,3,4}},{{0,4,5}}},
             1.0);

    QVERIFY(generate(input,createOptions(false),output));
    verifyMesh(input,output);
    QVERIFY(std::fabs(findTotalVolume(output) - 3.0) < 1.0e-10);
    QVERIFY(std::fabs(findTotalArea(output) - 14.0) < 1.0e-10);

    RTetraMeshOptions options = createOptions(true);
    options.setVolumeConstraint(0.01);
    QVERIFY(generate(input,options,output));
    verifyMesh(input,output);
    QVERIFY(std::fabs(findTotalVolume(output) - 3.0) < 1.0e-10);
    QVERIFY(std::fabs(findTotalArea(output) - 14.0) < 1.0e-10);
    QVERIFY(std::fabs(findTotalArea(output,2) - 8.0) < 1.0e-10);
    for (int i=0;i<int(output.getTetrahedra().size());i++)
    {
        QVERIFY(output.findTetrahedronVolume(i) <= 0.01 * (1.0 + 1.0e-12));
    }
}

void TestTetraMesher::sharpWedge()
{
    // Prism with a sharp edge, angle between its faces is 14 degrees.
    RTetraMeshData input, output;
    addPrism(input,
             {{{0.0,0.0}},{{4.0,0.0}},{{4.0,1.0}}},
             {{{0,1,2}}},
             1.0);

    RTetraMeshOptions options = createOptions(true);
    options.setVolumeConstraint(0.01);
    QVERIFY(generate(input,options,output));

    verifyMesh(input,output);
    QVERIFY(std::fabs(findTotalVolume(output) - 2.0) < 1.0e-10);
    for (int i=0;i<int(output.getTetrahedra().size());i++)
    {
        QVERIFY(output.findTetrahedronVolume(i) <= 0.01 * (1.0 + 1.0e-12));
    }
}

void TestTetraMesher::octahedron()
{
    // Shape which faces are not aligned with axes and do not form flat facets.
    RTetraMeshData input, output;
    input.addNode(1.0,0.0,0.0);
    input.addNode(-1.0,0.0,0.0);
    input.addNode(0.0,1.0,0.0);
    input.addNode(0.0,-1.0,0.0);
    input.addNode(0.0,0.0,1.0);
    input.addNode(0.0,0.0,-1.0);
    const int triangles[8][3] = {{0,2,4},{2,1,4},{1,3,4},{3,0,4},{2,0,5},{1,2,5},{3,1,5},{0,3,5}};
    for (const auto &triangle : triangles)
    {
        input.addTriangle(triangle[0],triangle[1],triangle[2],0);
    }

    RTetraMeshOptions options = createOptions(true);
    options.setVolumeConstraint(0.002);
    QVERIFY(generate(input,options,output));

    verifyMesh(input,output);
    QVERIFY(std::fabs(findTotalVolume(output) - 4.0/3.0) < 1.0e-10);
    QVERIFY(std::fabs(findTotalArea(output) - 4.0*std::sqrt(3.0)) < 1.0e-10);
}

void TestTetraMesher::sphere()
{
    // Curved surface, where no triangles lie in one plane.
    RTetraMeshData input, output;
    addSphere(input,2);

    // Volume enclosed by the surface.
    double volume = 0.0;
    for (const RTetraMeshData::Triangle &rTriangle : input.getTriangles())
    {
        const RTetraMeshData::Node &a = input.getNodes()[size_t(rTriangle.nodeIDs[0])];
        const RTetraMeshData::Node &b = input.getNodes()[size_t(rTriangle.nodeIDs[1])];
        const RTetraMeshData::Node &c = input.getNodes()[size_t(rTriangle.nodeIDs[2])];
        volume += (a.x*(b.y*c.z - b.z*c.y) + a.y*(b.z*c.x - b.x*c.z) + a.z*(b.x*c.y - b.y*c.x)) / 6.0;
    }
    const double area = findTotalArea(input);

    const double maxVolume = 0.005;
    const double maxRatio = 2.0;
    RTetraMeshOptions options = createOptions(true);
    options.setVolumeConstraint(maxVolume);
    options.setRadiusEdgeRatio(maxRatio);
    QVERIFY(generate(input,options,output));

    verifyMesh(input,output);
    QVERIFY(std::fabs(findTotalVolume(output) - std::fabs(volume)) < 1.0e-9);
    QVERIFY(std::fabs(findTotalArea(output) - area) < 1.0e-9);
    for (int i=0;i<int(output.getTetrahedra().size());i++)
    {
        QVERIFY(output.findTetrahedronVolume(i) <= maxVolume * (1.0 + 1.0e-12));
        QVERIFY(findRadiusEdgeRatio(output,i) <= maxRatio * (1.0 + 1.0e-12));
    }
}

void TestTetraMesher::sphereMinDihedralAngle()
{
    RTetraMeshData input, plain, improved;
    addSphere(input,3);

    RTetraMeshOptions options = createOptions(true);
    options.setVolumeConstraint(0.005);
    QVERIFY(generate(input,options,plain));

    const double minAngle = 10.0;
    options.setMinDihedralAngle(minAngle);
    QVERIFY(generate(input,options,improved));
    verifyMesh(input,improved);

    // Flat tetrahedra are found in a mesh refined only by radius/edge ratio, most of them have to be gone.
    int nPlainFlat = 0;
    for (int i=0;i<int(plain.getTetrahedra().size());i++)
    {
        if (findMinDihedralAngle(plain,i) < minAngle)
        {
            nPlainFlat++;
        }
    }
    int nImprovedFlat = 0;
    for (int i=0;i<int(improved.getTetrahedra().size());i++)
    {
        if (findMinDihedralAngle(improved,i) < minAngle)
        {
            nImprovedFlat++;
        }
    }
    // Optimization may already remove all flat cells without an explicit angle target.
    QVERIFY(nImprovedFlat * 10 <= nPlainFlat);
}

void TestTetraMesher::twoRegions()
{
    // Box divided by a wall into two volumes.
    RTetraMeshData input, output;
    addDividedBox(input,true);

    RTetraMeshOptions options = createOptions(true);
    options.setVolumeConstraint(0.01);
    QVERIFY(generate(input,options,output));

    verifyMesh(input,output);

    std::map<int,double> volumes;
    for (int i=0;i<int(output.getTetrahedra().size());i++)
    {
        volumes[output.getTetrahedra()[size_t(i)].marker] += output.findTetrahedronVolume(i);
    }
    QCOMPARE(int(volumes.size()),2);
    QVERIFY(std::fabs(volumes[0] - 1.0) < 1.0e-10);
    QVERIFY(std::fabs(volumes[1] - 1.0) < 1.0e-10);
    QVERIFY(std::fabs(findTotalArea(output,1) - 1.0) < 1.0e-10);
    QVERIFY(std::fabs(findTotalArea(output,0) - 10.0) < 1.0e-10);
}

void TestTetraMesher::nestedBoxes()
{
    // Box inside of a box is a volume of its own.
    RTetraMeshData input, output;
    addBox(input,0.0,0.0,0.0,3.0,3.0,3.0,0);
    addBox(input,1.0,1.0,1.0,2.0,2.0,2.0,1);

    RTetraMeshOptions options = createOptions(true);
    options.setVolumeConstraint(0.05);
    QVERIFY(generate(input,options,output));

    verifyMesh(input,output);

    std::map<int,double> volumes;
    for (int i=0;i<int(output.getTetrahedra().size());i++)
    {
        volumes[output.getTetrahedra()[size_t(i)].marker] += output.findTetrahedronVolume(i);
    }
    QCOMPARE(int(volumes.size()),2);
    std::vector<double> sorted;
    for (const auto &rVolume : volumes)
    {
        sorted.push_back(rVolume.second);
    }
    std::sort(sorted.begin(),sorted.end());
    QVERIFY(std::fabs(sorted[0] - 1.0) < 1.0e-9);
    QVERIFY(std::fabs(sorted[1] - 26.0) < 1.0e-9);
    QVERIFY(std::fabs(findTotalArea(output,1) - 6.0) < 1.0e-9);
}

void TestTetraMesher::pointsAndEdges()
{
    RTetraMeshData input, output;
    addBox(input,0.0,0.0,0.0,1.0,1.0,1.0,0);

    // Point on a corner and a point inside.
    input.getNodes()[0].marker = 4;
    const int inside = input.addNode(0.5,0.5,0.25,7);
    // Edge on an edge of the box and an edge inside.
    input.addEdge(0,1,2);
    const int n1 = input.addNode(0.25,0.5,0.5);
    const int n2 = input.addNode(0.75,0.5,0.5);
    input.addEdge(n1,n2,5);

    RTetraMeshOptions options = createOptions(true);
    options.setVolumeConstraint(0.002);
    QVERIFY(generate(input,options,output));

    verifyMesh(input,output);
    QVERIFY(std::fabs(findTotalVolume(output) - 1.0) < 1.0e-10);

    // Markers of points are not given to added nodes.
    QCOMPARE(output.getNodes()[0].marker,4);
    QCOMPARE(output.getNodes()[size_t(inside)].marker,7);
    int nMarked = 0;
    for (const RTetraMeshData::Node &rNode : output.getNodes())
    {
        if (rNode.marker != RTetraMeshData::noMarker)
        {
            nMarked++;
        }
    }
    QCOMPARE(nMarked,2);

    // Edges are split, but they keep their length and place.
    std::map<int,double> lengths;
    std::set<std::array<int,2>> tetrahedronEdges;
    for (const RTetraMeshData::Tetrahedron &rTetrahedron : output.getTetrahedra())
    {
        for (int i=0;i<4;i++)
        {
            for (int j=0;j<4;j++)
            {
                tetrahedronEdges.insert({{rTetrahedron.nodeIDs[i],rTetrahedron.nodeIDs[j]}});
            }
        }
    }
    for (const RTetraMeshData::Edge &rEdge : output.getEdges())
    {
        lengths[rEdge.marker] += findDistance(output,rEdge.nodeIDs[0],rEdge.nodeIDs[1]);
        QVERIFY(tetrahedronEdges.find({{rEdge.nodeIDs[0],rEdge.nodeIDs[1]}}) != tetrahedronEdges.end());
        for (int nodeID : rEdge.nodeIDs)
        {
            const RTetraMeshData::Node &rNode = output.getNodes()[size_t(nodeID)];
            if (rEdge.marker == 2)
            {
                QVERIFY(std::fabs(rNode.y) < 1.0e-12 && std::fabs(rNode.z) < 1.0e-12);
            }
            else
            {
                QVERIFY(std::fabs(rNode.y - 0.5) < 1.0e-12 && std::fabs(rNode.z - 0.5) < 1.0e-12);
            }
        }
    }
    QCOMPARE(int(lengths.size()),2);
    QVERIFY(std::fabs(lengths[2] - 1.0) < 1.0e-12);
    QVERIFY(std::fabs(lengths[5] - 0.5) < 1.0e-12);
}

void TestTetraMesher::quadrilateralFaces()
{
    // Sides of a box triangulated in a way which is not Delaunay.
    RTetraMeshData input, output;
    addBox(input,0.0,0.0,0.0,1.0,1.0,1.0,0);
    const int center = input.addNode(0.9,0.9,1.0);
    // Top side is replaced by a fan around a node close to its corner.
    std::vector<RTetraMeshData::Triangle> &rTriangles = input.getTriangles();
    rTriangles.erase(rTriangles.begin()+2,rTriangles.begin()+4);
    input.addTriangle(4,5,center,0);
    input.addTriangle(5,6,center,0);
    input.addTriangle(6,7,center,0);
    input.addTriangle(7,4,center,0);

    QVERIFY(generate(input,createOptions(false),output));

    verifyMesh(input,output);
    QVERIFY(std::fabs(findTotalVolume(output) - 1.0) < 1.0e-12);
    QVERIFY(std::fabs(findTotalArea(output) - 6.0) < 1.0e-12);
}

void TestTetraMesher::reconstruct()
{
    // Mesh of two volumes which are not separated by triangles.
    RTetraMeshData input, first, output;
    addDividedBox(input,true);

    RTetraMeshOptions options = createOptions(true);
    options.setVolumeConstraint(0.05);
    QVERIFY(generate(input,options,first));

    RTetraMeshData second;
    second.getNodes() = first.getNodes();
    for (const RTetraMeshData::Triangle &rTriangle : first.getTriangles())
    {
        if (rTriangle.marker == 0)
        {
            second.getTriangles().push_back(rTriangle);
        }
    }
    for (RTetraMeshData::Tetrahedron tetrahedron : first.getTetrahedra())
    {
        // Volume on the left side gets marker 5, the one on the right side gets 9.
        const RTetraMeshData::Node &rNode = first.getNodes()[size_t(tetrahedron.nodeIDs[0])];
        const RTetraMeshData::Node &rNode2 = first.getNodes()[size_t(tetrahedron.nodeIDs[1])];
        const RTetraMeshData::Node &rNode3 = first.getNodes()[size_t(tetrahedron.nodeIDs[2])];
        const RTetraMeshData::Node &rNode4 = first.getNodes()[size_t(tetrahedron.nodeIDs[3])];
        tetrahedron.marker = (rNode.x + rNode2.x + rNode3.x + rNode4.x) / 4.0 < 1.0 ? 5 : 9;
        second.getTetrahedra().push_back(tetrahedron);
    }

    options.setReconstruct(true);
    options.setVolumeConstraint(0.005);
    QVERIFY(generate(second,options,output));

    verifyMesh(second,output);
    QVERIFY(output.getTetrahedra().size() > 2 * first.getTetrahedra().size());

    std::map<int,double> volumes;
    for (int i=0;i<int(output.getTetrahedra().size());i++)
    {
        volumes[output.getTetrahedra()[size_t(i)].marker] += output.findTetrahedronVolume(i);
    }
    QCOMPARE(int(volumes.size()),2);
    QVERIFY(std::fabs(volumes[5] - 1.0) < 1.0e-10);
    QVERIFY(std::fabs(volumes[9] - 1.0) < 1.0e-10);
    // Interface between volumes has no marker.
    QVERIFY(std::fabs(findTotalArea(output,0) - 10.0) < 1.0e-10);
    QVERIFY(std::fabs(findTotalArea(output,RTetraMeshData::noMarker) - 1.0) < 1.0e-10);
}

void TestTetraMesher::duplicateNodes()
{
    // Two boxes which touch each other and do not share nodes.
    RTetraMeshData input, output;
    addBox(input,0.0,0.0,0.0,1.0,1.0,1.0,0);
    addBox(input,1.0,0.0,0.0,2.0,1.0,1.0,1,true);

    QVERIFY(generate(input,createOptions(false),output));

    QCOMPARE(output.getNodes().size(),input.getNodes().size());
    QVERIFY(std::fabs(findTotalVolume(output) - 2.0) < 1.0e-12);
    std::map<int,double> volumes;
    for (int i=0;i<int(output.getTetrahedra().size());i++)
    {
        volumes[output.getTetrahedra()[size_t(i)].marker] += output.findTetrahedronVolume(i);
    }
    QCOMPARE(int(volumes.size()),2);
}

void TestTetraMesher::invalidInput()
{
    RTetraMeshData output;

    // No triangles.
    {
        RTetraMeshData input;
        input.addNode(0.0,0.0,0.0);
        input.addNode(1.0,0.0,0.0);
        input.addNode(0.0,1.0,0.0);
        input.addNode(0.0,0.0,1.0);
        QVERIFY_THROWS_EXCEPTION(RError,RTetraMesher::generate(input,createOptions(false),output));
    }
    // Node which does not exist.
    {
        RTetraMeshData input;
        addBox(input,0.0,0.0,0.0,1.0,1.0,1.0,0);
        input.addTriangle(0,1,100,0);
        QVERIFY_THROWS_EXCEPTION(RError,RTetraMesher::generate(input,createOptions(false),output));
    }
    // Surface which is not closed.
    {
        RTetraMeshData input;
        addBox(input,0.0,0.0,0.0,1.0,1.0,1.0,0);
        input.getTriangles().pop_back();
        QVERIFY_THROWS_EXCEPTION(RError,RTetraMesher::generate(input,createOptions(false),output));
    }
}

QTEST_APPLESS_MAIN(TestTetraMesher)
#include "tst_tetra_mesher.moc"
