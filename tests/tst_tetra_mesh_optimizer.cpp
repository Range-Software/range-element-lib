#include <QTest>

#include <rbl_error.h>

#include <QTemporaryDir>
#include <QFile>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>
#include <set>

#include "rel_delaunay_tetrahedralizer.h"
#include "rel_geometric_predicates.h"
#include "rel_tetra_mesh_optimizer.h"
#include "rel_tetra_mesh_quality.h"

class TestTetraMeshOptimizer : public QObject
{
    Q_OBJECT
    static std::vector<int> finite(RTetraMesh &mesh)
    {
        std::vector<int> ids;
        for (int i = 0; i < mesh.getNTetraSlots(); i++)
        {
            if (mesh.isAlive(i) && !mesh.isGhost(i))
            {
                ids.push_back(i);
            }
        }
        return ids;
    }

    static void orient(RTetraMesh &mesh, std::vector<std::array<int, 4>> &tets)
    {
        for (auto &t : tets)
        {
            if (RGeometricPredicates::orient3d(mesh.getPoint(t[0]), mesh.getPoint(t[1]), mesh.getPoint(t[2]), mesh.getPoint(t[3])) < 0)
            {
                std::swap(t[0], t[1]);
            }
        }
    }

    static void bipyramid(RTetraMesh &mesh, double height = 1.0, double depth = -.1)
    {
        mesh.addPoint(0, 0, 0);
        mesh.addPoint(1, 0, 0);
        mesh.addPoint(0, 1, 0);
        mesh.addPoint(.2, .2, depth);
        mesh.addPoint(.2, .2, height);
        std::vector<int> duplicates;
        RDelaunayTetrahedralizer::tetrahedralize(mesh, duplicates);
        for (int id : finite(mesh))
        {
            mesh.setRegion(id, 0);
        }
    }

    static void hull(RBoundaryRecovery &boundary)
    {
        for (const auto &f : std::vector<std::array<int, 3>>{{0, 1, 3}, {1, 2, 3}, {2, 0, 3}, {0, 1, 4}, {1, 2, 4}, {2, 0, 4}})
        {
            boundary.addTriangle(f[0], f[1], f[2], 0);
        }
    }

    static void topology(const RTetraMesh &mesh)
    {
        for (int i = 0; i < mesh.getNTetraSlots(); i++)
        {
            if (mesh.isAlive(i))
            {
                const auto &t = mesh.getTetra(i);
                if (!mesh.isGhost(i))
                {
                    QVERIFY(RGeometricPredicates::orient3d(mesh.getPoint(t.v[0]), mesh.getPoint(t.v[1]), mesh.getPoint(t.v[2]), mesh.getPoint(t.v[3])) > 0);
                }
                for (int j = 0; j < 4; j++)
                {
                    QVERIFY(t.n[j] >= 0 && t.n[j] < mesh.getNTetraSlots());
                    QVERIFY(mesh.isAlive(t.n[j]));
                    const auto &n = mesh.getTetra(t.n[j]);
                    bool back = false;
                    for (int k = 0; k < 4; k++)
                    {
                        if (n.n[k] == i)
                        {
                            std::array<int, 3> a, b;
                            for (int x = 0; x < 3; x++)
                            {
                                a[x] = t.v[RTetraMesh::faceVertices[j][x]];
                                b[x] = n.v[RTetraMesh::faceVertices[k][x]];
                            }
                            std::sort(a.begin(), a.end());
                            std::sort(b.begin(), b.end());
                            QCOMPARE(a, b);
                            back = true;
                        }
                    }
                    QVERIFY(back);
                }
            }
        }
        for (int p = 0; p < mesh.getNPoints(); p++)
        {
            int id = mesh.getVertexTetra(p);
            QVERIFY(id >= 0 && mesh.isAlive(id));
            const auto &t = mesh.getTetra(id);
            QVERIFY(std::find(t.v, t.v + 4, p) != t.v + 4);
        }
    }

    private slots:

        void planarPatchRecovery_data()
        {
            QTest::addColumn<double>("scale");
            QTest::addColumn<bool>("protectDiagonal");
            QTest::newRow("small") << 1e-6 << false;
            QTest::newRow("unit") << 1.0 << false;
            QTest::newRow("large") << 1e6 << false;
            QTest::newRow("constrained-diagonal") << 1.0 << true;
        }

        void planarPatchRecovery()
        {
            QFETCH(double, scale);
            QFETCH(bool, protectDiagonal);
            RTetraMesh mesh;
            for (const auto &p : std::vector<std::array<double, 3>>{
                     {0, 0, 0}, {3, 0, 0}, {4, 2, 0}, {2, 4, 0}, {0, 2, 0}, {.8, 1.3, 1}, {.8, 1.3, -1}})
            {
                mesh.addPoint(p[0] * scale, p[1] * scale, p[2] * scale);
            }
            std::vector<int> duplicates;
            RDelaunayTetrahedralizer::tetrahedralize(mesh, duplicates);
            auto fan = [&](int root)
            {
                std::vector<std::array<int, 4>> cells;
                for (int i = 1; i < 4; i++)
                {
                    for (int pole : {5, 6})
                    {
                        cells.push_back({root, (root + i) % 5, (root + i + 1) % 5, pole});
                    }
                }
                orient(mesh, cells);
                return cells;
            };
            auto initial = fan(0);
            QVERIFY(mesh.replaceTetras(finite(mesh), initial));
            RBoundaryRecovery boundary(mesh, 1e-10);
            for (int i = 0; i < 5; i++)
            {
                for (int pole : {5, 6})
                {
                    boundary.addTriangle(i, (i + 1) % 5, pole, 1);
                }
            }
            for (int i = 1; i < 4; i++)
            {
                boundary.addTriangle(0, i, i + 1, 7);
            }
            if (protectDiagonal)
            {
                boundary.addEdge(0, 2, 42);
            }
            boundary.initialize();
            // Replace a recovered planar fan by a different, non-cocircular fan.
            // A single diagonal flip cannot match these two triangulations.
            auto replacement = fan(1);
            QVERIFY(mesh.replaceTetras(finite(mesh), replacement));
            boundary.setMinLength(scale * 1e-8);
            boundary.setMaxSteinerVertices(500);
            boundary.recover();
            if (!protectDiagonal)
            {
                QCOMPARE(mesh.getNPoints(), 7);
            }
            double interfaceArea = 0;
            for (const auto &face : boundary.getSubFaces())
            {
                if (face.alive)
                {
                    RTetraMesh::Face found;
                    QVERIFY(mesh.findFace(face.v[0], face.v[1], face.v[2], found));
                    if (face.marker == 7)
                    {
                        const double *a = mesh.getPoint(face.v[0]), *b = mesh.getPoint(face.v[1]), *c = mesh.getPoint(face.v[2]);
                        QCOMPARE(a[2], 0.0);
                        QCOMPARE(b[2], 0.0);
                        QCOMPARE(c[2], 0.0);
                        interfaceArea += std::fabs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) / 2;
                    }
                    else
                    {
                        QCOMPARE(face.marker, 1);
                    }
                }
            }
            QVERIFY(std::fabs(interfaceArea / (scale * scale) - 11.0) < 1e-10);
            if (protectDiagonal)
            {
                double length = 0;
                for (const auto &edge : boundary.getSubSegments())
                {
                    if (edge.alive && edge.marker == 42)
                    {
                        QVERIFY(mesh.findEdge(edge.v[0], edge.v[1]) != RTetraMesh::noTetra);
                        const double *a = mesh.getPoint(edge.v[0]), *b = mesh.getPoint(edge.v[1]);
                        length += std::hypot(a[0] - b[0], a[1] - b[1]);
                    }
                }
                QVERIFY(std::fabs(length - std::sqrt(20.0)) < 1e-10);
            }
            topology(mesh);
        }

        void recoveryFailureReportsSplitLimit()
        {
            RTetraMesh mesh;
            bipyramid(mesh);
            std::vector<std::array<int, 4>> cells{{3, 4, 0, 1}, {3, 4, 1, 2}, {3, 4, 2, 0}};
            orient(mesh, cells);
            QVERIFY(mesh.replaceTetras(finite(mesh), cells));
            RBoundaryRecovery boundary(mesh, 1e-10);
            hull(boundary);
            boundary.addTriangle(0, 1, 2, 8);
            boundary.initialize();
            boundary.setMinLength(10.0);
            const int points = mesh.getNPoints();
            bool failed = false;
            try
            {
                boundary.recover();
            }
            catch (const RError &error)
            {
                failed = true;
                const QString message = error.getMessage();
                QVERIFY2(message.contains("Triangle circumradius is at or below the minimum split length"), qPrintable(message));
                QVERIFY(message.contains("nodes 0, 1, 2"));
                QVERIFY(message.contains("surface marker 8"));
                QVERIFY(message.contains("Circumradius ="));
                QVERIFY(message.contains("minimum split length = 10"));
                QVERIFY(!message.contains("Surface may intersect itself"));
            }
            QVERIFY(failed);
            QCOMPARE(mesh.getNPoints(), points);
            topology(mesh);
        }

        void qualityAndScale()
        {
            for (double scale : {1e-8, 1.0, 1e8})
            {
                const double a[3] = {0, 0, 0}, b[3] = {scale, 0, 0}, c[3] = {.5 * scale, std::sqrt(3.0) * .5 * scale, 0}, d[3] = {.5 * scale, std::sqrt(3.0) / 6 * scale, std::sqrt(2.0 / 3) * scale};
                const auto q = RTetraMeshQuality::measure(a, b, c, d);
                QVERIFY(std::fabs(q.meanRatio - 1) < 1e-12);
                QVERIFY(std::fabs(q.radiusEdge - std::sqrt(6.0) / 4) < 1e-12);
                QVERIFY(std::fabs(q.minAngle - 70.5287793655) < 1e-8);
                double flat[3] = {.3 * scale, .2 * scale, 1e-15 * scale};
                const auto f = RTetraMeshQuality::measure(a, b, c, flat);
                QVERIFY(f.volume > 0);
                QVERIFY(f.scaledVolume < 1e-12);
                QVERIFY(f.meanRatio < 1e-8);
            }
        }

        void rejectedReplacementIsAtomic()
        {
            RTetraMesh mesh;
            bipyramid(mesh);
            const auto ids = finite(mesh);
            const uint64_t count = mesh.getNTetras();
            QVERIFY(!mesh.replaceTetras(ids, {{0, 1, 2, 3}}));
            QVERIFY(!mesh.replaceTetras({-1}, {{0, 1, 2, 3}}));
            QCOMPARE(mesh.getNTetras(), count);
            QCOMPARE(finite(mesh), ids);
            topology(mesh);
        }

        void replacementTrialRollback()
        {
            RTetraMesh mesh;
            bipyramid(mesh);
            const auto original = mesh;
            auto same = [&](const RTetraMesh &expected)
            {
                QCOMPARE(mesh.getNTetraSlots(), expected.getNTetraSlots());
                QCOMPARE(mesh.getNTetras(), expected.getNTetras());
                QCOMPARE(mesh.getCreatedTetras(), expected.getCreatedTetras());
                QCOMPARE(mesh.getCavityTetras(), expected.getCavityTetras());
                for (int i = 0; i < mesh.getNTetraSlots(); i++)
                {
                    const auto &a = mesh.getTetra(i), &b = expected.getTetra(i);
                    for (int j = 0; j < 4; j++)
                    {
                        QCOMPARE(a.v[j], b.v[j]);
                        QCOMPARE(a.n[j], b.n[j]);
                    }
                    QCOMPARE(a.region, b.region);
                    QCOMPARE(mesh.isInCavity(i), expected.isInCavity(i));
                }
                for (int i = 0; i < mesh.getNPoints(); i++)
                {
                    QCOMPARE(mesh.getVertexTetra(i), expected.getVertexTetra(i));
                    std::vector<int> a, b;
                    mesh.findStar(i, a);
                    expected.findStar(i, b);
                    QCOMPARE(a, b);
                }
                topology(mesh);
            };
            std::vector<std::array<int, 4>> three = {{3, 4, 0, 1}, {3, 4, 1, 2}, {3, 4, 2, 0}};
            orient(mesh, three);
            std::vector<std::array<int, 4>> two;
            for (int id : finite(mesh))
            {
                const auto &t = mesh.getTetra(id);
                two.push_back({t.v[0], t.v[1], t.v[2], t.v[3]});
            }
            // Appended slots, reused slots, multiple replacements, nested rollback,
            // and a rejected proposal all restore the exact original slot/hint order.
            for (int repeat = 0; repeat < 2; repeat++)
            {
                RTetraMesh::ReplacementTrial trial(mesh);
                QVERIFY(!trial.replace({-1}, three));
                same(original);
                QVERIFY(trial.replace(finite(mesh), three));
                topology(mesh);
                QVERIFY(trial.bytesSaved() > 0);
                const auto intermediate = mesh;
                {
                    RTetraMesh::ReplacementTrial inner(mesh);
                    QVERIFY(inner.replace(finite(mesh), two));
                    topology(mesh);
                }
                same(intermediate);
                QVERIFY(trial.replace(finite(mesh), two));
                topology(mesh);
                trial.rollback();
                same(original);
            }
            try
            {
                RTetraMesh::ReplacementTrial trial(mesh);
                QVERIFY(trial.replace(finite(mesh), three));
                throw 7;
            }
            catch (int)
            {
                same(original);
            }
            // A real replacement after rollback must allocate exactly as a fresh copy.
            auto expected = original;
            QVERIFY(expected.replaceTetras(finite(expected), three));
            QVERIFY(mesh.replaceTetras(finite(mesh), three));
            same(expected);
        }

        void optimizerRefreshesExternalChanges()
        {
            RTetraMesh mesh;
            bipyramid(mesh);
            RBoundaryRecovery boundary(mesh, 1e-10);
            hull(boundary);
            boundary.initialize();
            boundary.recover();
            boundary.classifyRegions();
            RTetraMeshOptions options;
            RTetraMeshOptimizer reused(mesh, boundary, options);
            reused.optimize(false);
            const double moved[3] = {.25, .18, 1.1};
            mesh.setPoint(4, moved);
            auto expected = mesh;
            RBoundaryRecovery expectedBoundary(expected, boundary);
            RTetraMeshOptimizer fresh(expected, expectedBoundary, options);
            fresh.optimize(false);
            reused.optimize(false);
            QCOMPARE(mesh.getNTetraSlots(), expected.getNTetraSlots());
            QCOMPARE(mesh.getNTetras(), expected.getNTetras());
            for (int i = 0; i < mesh.getNTetraSlots(); i++)
            {
                const auto &a = mesh.getTetra(i), &b = expected.getTetra(i);
                for (int j = 0; j < 4; j++)
                {
                    QCOMPARE(a.v[j], b.v[j]);
                    QCOMPARE(a.n[j], b.n[j]);
                }
                QCOMPARE(a.region, b.region);
            }
            topology(mesh);
        }

        void parallelQualitySummary()
        {
            RTetraMeshData data;
            data.addNode(0, 0, 0);
            data.addNode(1, 0, 0);
            data.addNode(0, 1, 0);
            data.addNode(.2, .2, .001);
            data.addTriangle(0, 1, 2, 7);
            data.addTriangle(0, 1, 3, 7);
            RTetraMeshData::Tetrahedron cell;
            std::copy_n(std::array<int, 4>{0, 1, 2, 3}.begin(), 4, cell.nodeIDs);
            // Match the library's positive orientation convention.
            const double a[3] = {0, 0, 0}, b[3] = {1, 0, 0}, c[3] = {0, 1, 0}, d[3] = {.2, .2, .001};
            if (RGeometricPredicates::orient3d(a, b, c, d) < 0)
            {
                std::swap(cell.nodeIDs[0], cell.nodeIDs[1]);
            }
            for (int i = 0; i < 32769; i++)
            {
                auto t = cell;
                if (i % 31 == 0)
                {
                    t.nodeIDs[0] = -1;
                }
                data.getTetrahedra().push_back(t);
            }
            RTetraMeshOptions options;
            options.setUseSizeFunction(true);
            options.setMinDihedralAngle(10);
            options.setVolumeConstraint(1e-5);
            for (auto &v : data.getNodes())
            {
                v.size = .01;
            }
            const auto serial = RTetraMeshQuality::summarize(data, options, 1);
            for (unsigned workers : {0u, 2u, 4u})
            {
                const auto parallel = RTetraMeshQuality::summarize(data, options, workers);
                QCOMPARE(parallel.elements, serial.elements);
                QCOMPARE(parallel.invalid, serial.invalid);
                QCOMPARE(parallel.nearFlat, serial.nearFlat);
                QCOMPARE(parallel.badRadius, serial.badRadius);
                QCOMPARE(parallel.badAngle, serial.badAngle);
                QCOMPARE(parallel.constrainedAngle, serial.constrainedAngle);
                QCOMPARE(parallel.tooLarge, serial.tooLarge);
                QCOMPARE(parallel.freeBelowOne, serial.freeBelowOne);
                QCOMPARE(parallel.prescribedBelowOne, serial.prescribedBelowOne);
                QCOMPARE(parallel.minAngle, serial.minAngle);
                QCOMPARE(parallel.radiusEdge95, serial.radiusEdge95);
                QCOMPARE(parallel.maxRadiusEdge, serial.maxRadiusEdge);
            }
            QCOMPARE(RTetraMeshQuality::summarize(RTetraMeshData(), options, 4).elements, size_t(0));
        }

        void constrainedFlip_data()
        {
            QTest::addColumn<int>("kind");
            QTest::newRow("2-to-3") << 0;
            QTest::newRow("protected-face") << 1;
            QTest::newRow("region-interface") << 2;
            QTest::newRow("3-to-2") << 3;
            QTest::newRow("protected-edge") << 4;
        }

        void constrainedFlip()
        {
            QFETCH(int, kind);
            const bool three = kind >= 3;
            RTetraMesh mesh;
            bipyramid(mesh, 1.0, three ? -1.0 : -.1);
            std::vector<std::array<int, 4>> bad = three
                                                      ? std::vector<std::array<int, 4>>{{3, 4, 0, 1}, {3, 4, 1, 2}, {3, 4, 2, 0}}
                                                      : std::vector<std::array<int, 4>>{{0, 1, 2, 3}, {0, 2, 1, 4}};
            orient(mesh, bad);
            QVERIFY(mesh.replaceTetras(finite(mesh), bad));
            RBoundaryRecovery boundary(mesh, 1e-10);
            hull(boundary);
            if (kind == 1)
            {
                boundary.addTriangle(0, 1, 2, 1);
            }
            if (kind == 4)
            {
                boundary.addEdge(3, 4, 7);
            }
            boundary.initialize();
            boundary.recover();
            auto ids = finite(mesh);
            if (kind == 2)
            {
                mesh.setRegion(ids.front(), 1);
            }
            boundary.getCreatedTetras() = ids;
            RTetraMeshOptions options;
            RTetraMeshOptimizer optimizer(mesh, boundary, options);
            const auto changes = optimizer.repairCreated();
            if (kind == 1 || kind == 2 || kind == 4)
            {
                QCOMPARE(changes, 0u);
                QCOMPARE(finite(mesh).size(), bad.size());
            }
            else
            {
                QVERIFY(changes > 0);
                QCOMPARE(finite(mesh).size(), size_t(three ? 2 : 3));
            }
            for (const auto &face : boundary.getSubFaces())
            {
                if (face.alive)
                {
                    RTetraMesh::Face found;
                    QVERIFY(mesh.findFace(face.v[0], face.v[1], face.v[2], found));
                }
            }
            topology(mesh);
        }

        void fourToFour_data()
        {
            QTest::addColumn<bool>("protect");
            QTest::newRow("repair") << false;
            QTest::newRow("protected") << true;
        }

        void fourToFour()
        {
            QFETCH(bool, protect);
            RTetraMesh mesh;
            for (const auto &p : std::vector<std::array<double, 3>>{{.1, 0, -1}, {.1, 0, 1}, {-1, 0, 0}, {0, -.5, 0}, {1, 0, 0}, {0, .7, 0}})
            {
                mesh.addPoint(p[0], p[1], p[2]);
            }
            std::vector<int> duplicates;
            RDelaunayTetrahedralizer::tetrahedralize(mesh, duplicates);
            for (int id : finite(mesh))
            {
                mesh.setRegion(id, 0);
            }
            std::vector<std::array<int, 4>> bad = {{0, 1, 2, 3}, {0, 1, 3, 4}, {0, 1, 4, 5}, {0, 1, 5, 2}};
            orient(mesh, bad);
            QVERIFY(mesh.replaceTetras(finite(mesh), bad));
            RBoundaryRecovery boundary(mesh, 1e-10);
            for (int i = 0; i < 4; i++)
            {
                for (int pole = 0; pole < 2; pole++)
                {
                    boundary.addTriangle(pole, 2 + i, 2 + (i + 1) % 4, 0);
                }
            }
            if (protect)
            {
                boundary.addEdge(0, 1, 7);
            }
            boundary.initialize();
            boundary.recover();
            boundary.getCreatedTetras() = finite(mesh);
            RTetraMeshOptions options;
            RTetraMeshOptimizer optimizer(mesh, boundary, options);
            const auto changes = optimizer.repairCreated();
            if (protect)
            {
                QVERIFY(mesh.findEdge(0, 1) >= 0);
            }
            else
            {
                QVERIFY(changes > 0);
                QCOMPARE(mesh.findEdge(0, 1), RTetraMesh::noTetra);
                QCOMPARE(finite(mesh).size(), size_t(4));
            }
            for (const auto &face : boundary.getSubFaces())
            {
                if (face.alive)
                {
                    RTetraMesh::Face found;
                    QVERIFY(mesh.findFace(face.v[0], face.v[1], face.v[2], found));
                }
            }
            topology(mesh);
        }

        void nonFanLink()
        {
            RTetraMesh mesh;
            mesh.addPoint(0, 0, -1.3);
            mesh.addPoint(0, 0, 1.3);
            for (int i = 0; i < 6; i++)
            {
                mesh.addPoint(std::cos(i * std::acos(-1.0) / 3), std::sin(i * std::acos(-1.0) / 3), 0);
            }
            std::vector<int> duplicates;
            RDelaunayTetrahedralizer::tetrahedralize(mesh, duplicates);
            for (int id : finite(mesh))
            {
                mesh.setRegion(id, 0);
            }
            std::vector<std::array<int, 4>> bad;
            for (int i = 0; i < 6; i++)
            {
                bad.push_back({0, 1, 2 + i, 2 + (i + 1) % 6});
            }
            orient(mesh, bad);
            QVERIFY(mesh.replaceTetras(finite(mesh), bad));
            RBoundaryRecovery boundary(mesh, 1e-10);
            for (int i = 0; i < 6; i++)
            {
                for (int pole = 0; pole < 2; pole++)
                {
                    boundary.addTriangle(pole, 2 + i, 2 + (i + 1) % 6, 0);
                }
            }
            boundary.initialize();
            boundary.recover();
            boundary.getCreatedTetras() = finite(mesh);
            RTetraMeshOptions options;
            RTetraMeshOptimizer optimizer(mesh, boundary, options);
            QVERIFY(optimizer.repairCreated() > 0);
            QCOMPARE(mesh.findEdge(0, 1), RTetraMesh::noTetra);
            // No link vertex occurs in every new link triangle: this cannot be a fan.
            std::set<int> common = {2, 3, 4, 5, 6, 7};
            for (int id : finite(mesh))
            {
                const auto &t = mesh.getTetra(id);
                for (auto it = common.begin(); it != common.end();)
                {
                    if (std::find(t.v, t.v + 4, *it) == t.v + 4)
                    {
                        it = common.erase(it);
                    }
                    else
                    {
                        ++it;
                    }
                }
            }
            QVERIFY(common.empty());
            QCOMPARE(finite(mesh).size(), size_t(8));
            topology(mesh);
        }

        void boundedOffCenters_data()
        {
            QTest::addColumn<unsigned int>("nodeBudget");
            QTest::newRow("bounded-default") << 0u;
            QTest::newRow("user-limit") << 3u;
        }

        void boundedOffCenters()
        {
            QFETCH(unsigned int, nodeBudget);
            RTetraMesh mesh;
            for (int i = 0; i < 8; i++)
            {
                mesh.addPoint(i & 1, (i >> 1) & 1, (i >> 2) & 1);
            }
            uint64_t state = 123;
            auto random = [&]()
            {
                state = state * 6364136223846793005ULL + 1442695040888963407ULL;
                return double(state >> 11) / double(1ULL << 53);
            };
            for (int i = 0; i < 96; i++)
            {
                const double x = random(), y = random(), z = random();
                mesh.addPoint(x, y, z);
            }
            const int before = mesh.getNPoints();
            std::vector<int> duplicates;
            RDelaunayTetrahedralizer::tetrahedralize(mesh, duplicates);
            RBoundaryRecovery boundary(mesh, 1e-10);
            for (int id : finite(mesh))
            {
                mesh.setRegion(id, 0);
                const auto &t = mesh.getTetra(id);
                for (int i = 0; i < 4; i++)
                {
                    if (mesh.isGhost(t.n[i]))
                    {
                        const int *f = RTetraMesh::faceVertices[i];
                        boundary.addTriangle(t.v[f[0]], t.v[f[1]], t.v[f[2]], 0);
                    }
                }
            }
            boundary.initialize();
            boundary.recover();
            RTetraMeshOptions options;
            options.setRadiusEdgeRatio(1.4);
            options.setMinDihedralAngle(10);
            options.setMaxSteinerNodes(nodeBudget);
            RTetraMeshOptimizer optimizer(mesh, boundary, options);
            optimizer.optimize();
            qInfo("off-centers %zu visits %zu", optimizer.getInsertedPoints(), optimizer.getQueueVisits());
            QVERIFY(optimizer.getInsertedPoints() > 0);
            QVERIFY(optimizer.getInsertedPoints() <= (nodeBudget > 0 ? nodeBudget : 16));
            QCOMPARE(mesh.getNPoints(), before + int(optimizer.getInsertedPoints() + optimizer.getBoundaryPoints()));
            if (nodeBudget > 0)
            {
                QVERIFY(optimizer.getInsertedPoints() + optimizer.getBoundaryPoints() <= nodeBudget);
            }
            QVERIFY(optimizer.getQueueVisits() < size_t(10000));
            for (const auto &face : boundary.getSubFaces())
            {
                if (face.alive)
                {
                    RTetraMesh::Face found;
                    QVERIFY(mesh.findFace(face.v[0], face.v[1], face.v[2], found));
                }
            }
            topology(mesh);
        }

        void patchReport()
        {
            RTetraMeshData data;
            data.addNode(0, 0, 0);
            data.addNode(1, 0, 0);
            data.addNode(.5, .01, 0);
            data.addNode(.5, .005, .001);
            data.addTetrahedron(0, 1, 2, 3);
            for (int i = 0; i < 4; i++)
            {
                const int *f = RTetraMesh::faceVertices[i];
                data.addTriangle(f[0], f[1], f[2]);
            }
            RTetraMeshQuality::PatchAttempts attempts;
            attempts.vertices = {0, 1, 2, 3};
            attempts.visits = 2;
            attempts.spacing = 3;
            attempts.boundary = 4;
            QTemporaryDir dir;
            QVERIFY(dir.isValid());
            const QString path = dir.filePath("patches.json");
            RTetraMeshOptions options;
            options.setMinDihedralAngle(10);
            RTetraMeshQuality::writeReport(data, options, {attempts}, path.toStdString());
            QFile file(path);
            QVERIFY(file.open(QIODevice::ReadOnly));
            const auto report = QJsonDocument::fromJson(file.readAll()).object();
            QCOMPARE(report["schemaVersion"].toInt(), 1);
            const auto patches = report["patches"].toArray();
            QCOMPARE(patches.size(), 1);
            const auto patch = patches[0].toObject();
            QCOMPARE(patch["coordinates"].toArray().size(), 4);
            QCOMPARE(patch["attempts"].toObject()["spacing"].toInt(), 3);
            QCOMPARE(patch["attempts"].toObject()["boundary"].toInt(), 4);
            QVERIFY(patch["minimumFreeDihedral"].toDouble() > patch["minimumDihedral"].toDouble());
            QVERIFY(patch["minimumOnPrescribedFaces"].toBool());
            QVERIFY(!patch["observations"].toArray().empty());
        }

        void independentBoundaryTrial()
        {
            RTetraMesh mesh;
            mesh.addPoint(0, 0, 0);
            mesh.addPoint(1, 0, 0);
            mesh.addPoint(0, 1, 0);
            mesh.addPoint(0, 0, 1);
            mesh.createFirstTetra(0, 1, 2, 3);
            RBoundaryRecovery boundary(mesh, 1e-10);
            for (int i = 0; i < 4; i++)
            {
                const int *f = RTetraMesh::faceVertices[i];
                boundary.addTriangle(f[0], f[1], f[2], 0);
            }
            boundary.initialize();
            boundary.recover();
            QCOMPARE(boundary.classifyRegions(), 1);
            RTetraMesh trial = mesh;
            RBoundaryRecovery trialBoundary(trial, boundary);
            QVERIFY(trialBoundary.splitSubSegment(0, .1));
            trialBoundary.recover();
            QCOMPARE(mesh.getNPoints(), 4);
            QCOMPARE(boundary.getNSteinerVertices(), uint64_t(0));
            QVERIFY(trial.getNPoints() > 4);
            QVERIFY(trialBoundary.getNSteinerVertices() > 0);
            QCOMPARE(trialBoundary.classifyRegions(), 1);
            QCOMPARE(boundary.classifyRegions(), 1);
            topology(mesh);
            topology(trial);
        }

        void insertionHistory()
        {
            RTetraMesh mesh;
            mesh.addPoint(0, 0, 0);
            mesh.addPoint(1, 0, 0);
            mesh.addPoint(0, 1, 0);
            mesh.addPoint(0, 0, 1);
            mesh.createFirstTetra(0, 1, 2, 3);
            RBoundaryRecovery boundary(mesh, 1e-10);
            const double first[3] = {.2, .2, .2};
            QVERIFY(boundary.prepareVolumeCavity(first, mesh.locate(first)));
            boundary.insertCavityPoint(first, RBoundaryRecovery::VertexVolume, .1);
            const double close[3] = {.201, .2, .2};
            QVERIFY(boundary.prepareVolumeCavity(close, mesh.locate(close)));
            // A smaller current-cell bound must not forget an earlier insertion scale.
            QVERIFY(!boundary.insertionAllowed(close, RBoundaryRecovery::VertexVolume, 1e-6));
            const double far[3] = {.4, .2, .2};
            QVERIFY(boundary.prepareVolumeCavity(far, mesh.locate(far)));
            QVERIFY(boundary.insertionAllowed(far, RBoundaryRecovery::VertexVolume, 1e-6));
            topology(mesh);
        }

        void smoothingPreservesBoundary_data()
        {
            QTest::addColumn<bool>("adopted");
            QTest::newRow("interior") << false;
            QTest::newRow("adopted-by-facet") << true;
        }

        void smoothingPreservesBoundary()
        {
            QFETCH(bool, adopted);
            RTetraMesh mesh;
            mesh.addPoint(0, 0, 0);
            mesh.addPoint(1, 0, 0);
            mesh.addPoint(0, 1, 0);
            mesh.addPoint(0, 0, 1);
            mesh.createFirstTetra(0, 1, 2, 3);
            RBoundaryRecovery boundary(mesh, 1e-10);
            for (int i = 0; i < 4; i++)
            {
                const int *f = RTetraMesh::faceVertices[i];
                boundary.addTriangle(f[0], f[1], f[2], 0);
            }
            boundary.initialize();
            boundary.recover();
            boundary.classifyRegions();
            const double p[3] = {.02, .03, .04};
            QVERIFY(boundary.prepareVolumeCavity(p, mesh.locate(p)));
            const int v = boundary.insertCavityPoint(p, RBoundaryRecovery::VertexVolume);
            boundary.recover();
            boundary.classifyRegions();
            if (adopted)
            {
                boundary.addTriangle(0, 1, v, 1);
            }
            RTetraMeshOptions options;
            RTetraMeshOptimizer optimizer(mesh, boundary, options);
            const auto changes = optimizer.optimize();
            if (adopted)
            {
                QVERIFY(boundary.isBoundaryVertex(v));
                for (int i = 0; i < 3; i++)
                {
                    QCOMPARE(mesh.getPoint(v)[i], p[i]);
                }
            }
            else
            {
                QVERIFY(changes > 0);
                QVERIFY(mesh.getPoint(v)[0] != p[0] || mesh.getPoint(v)[1] != p[1] || mesh.getPoint(v)[2] != p[2]);
            }
            const double fixed[4][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    QCOMPARE(mesh.getPoint(i)[j], fixed[i][j]);
                }
            }
            topology(mesh);
        }

        void diagnosticConstraints()
        {
            RTetraMeshData data;
            data.addNode(0, 0, 0);
            data.addNode(1, 0, 0);
            data.addNode(.5, .01, 0);
            data.addNode(.5, .005, .001);
            data.addTetrahedron(0, 1, 2, 3);
            RTetraMeshOptions options;
            options.setMinDihedralAngle(60);
            const auto free = RTetraMeshQuality::summarize(data, options);
            QCOMPARE(free.badAngle, size_t(1));
            QVERIFY(!free.targetsMet());
            for (int i = 0; i < 4; i++)
            {
                const int *f = RTetraMesh::faceVertices[i];
                data.addTriangle(f[0], f[1], f[2]);
            }
            const auto fixed = RTetraMeshQuality::summarize(data, options);
            QCOMPARE(fixed.badAngle, size_t(0));
            QCOMPARE(fixed.constrainedAngle, size_t(1));
            data.getTetrahedra()[0].nodeIDs[0] = -1;
            const auto invalid = RTetraMeshQuality::summarize(data, options);
            QCOMPARE(invalid.invalid, size_t(1));
            QVERIFY(!invalid.targetsMet());
        }
};
QTEST_APPLESS_MAIN(TestTetraMeshOptimizer)
#include "tst_tetra_mesh_optimizer.moc"
