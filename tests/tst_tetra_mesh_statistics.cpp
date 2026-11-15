#include <QTest>
#include <QScopeGuard>
#include <QRegularExpression>

#include <cmath>
#include <numeric>

#include <rbl_logger.h>

#include "rel_tetra_mesh.h"
#include "rel_tetra_mesh_statistics.h"

class TestTetraMeshStatistics : public QObject
{
    Q_OBJECT

    static RTetraMeshData regular(double scale = 1.0)
    {
        RTetraMeshData data;
        data.addNode(0, 0, 0);
        data.addNode(scale, 0, 0);
        data.addNode(scale / 2, std::sqrt(3.0) * scale / 2, 0);
        data.addNode(scale / 2, std::sqrt(3.0) * scale / 6, std::sqrt(2.0 / 3.0) * scale);
        data.addTetrahedron(0, 1, 2, 3, 0);
        for (int i = 0; i < 4; i++)
        {
            const int *face = RTetraMesh::faceVertices[i];
            data.addTriangle(face[0], face[1], face[2], 0);
        }
        return data;
    }

    template <size_t N>
    static uint64_t total(const std::array<uint64_t, N> &values)
    {
        return std::accumulate(values.begin(), values.end(), uint64_t(0));
    }

    private slots:

        void regularTetrahedron_data()
        {
            QTest::addColumn<double>("scale");
            QTest::newRow("small") << 1e-6;
            QTest::newRow("unit") << 1.0;
            QTest::newRow("large") << 1e6;
        }

        void regularTetrahedron()
        {
            QFETCH(double, scale);
            const auto data = regular(scale);
            const auto result = RTetraMeshStatistics::measure(data, data);
            QCOMPARE(result.validElements, size_t(1));
            QCOMPARE(result.invalidElements, size_t(0));
            QCOMPARE(result.faces, size_t(4));
            QCOMPARE(result.exteriorFaces, size_t(4));
            QVERIFY(std::fabs(result.volume.minimum / (scale * scale * scale) - std::sqrt(2.0) / 12) < 1e-12);
            QCOMPARE(result.volume.minimum, result.volume.maximum);
            QVERIFY(std::fabs(result.edge.minimum / scale - 1) < 1e-12);
            QVERIFY(std::fabs(result.edge.maximum / scale - 1) < 1e-12);
            QVERIFY(std::fabs(result.aspectRatio.minimum - std::sqrt(1.5)) < 1e-12);
            QVERIFY(std::fabs(result.aspectRatio.maximum - std::sqrt(1.5)) < 1e-12);
            QVERIFY(std::fabs(result.faceAngle.minimum - 60) < 1e-10);
            QVERIFY(std::fabs(result.faceAngle.maximum - 60) < 1e-10);
            QVERIFY(std::fabs(result.dihedralAngle.minimum - 70.5287793655093) < 1e-10);
            QVERIFY(std::fabs(result.dihedralAngle.maximum - 70.5287793655093) < 1e-10);
            QVERIFY(std::fabs(result.minimumInputFaceAngle - 60) < 1e-10);
            QVERIFY(std::fabs(result.minimumInputDihedralAngle - 70.5287793655093) < 1e-10);
            QCOMPARE(total(result.aspectHistogram), uint64_t(1));
            QCOMPARE(result.aspectHistogram[0], uint64_t(1));
            QCOMPARE(total(result.faceAngleHistogram), uint64_t(12));
            QCOMPARE(total(result.dihedralAngleHistogram), uint64_t(6));
            QCOMPARE(result.dihedralAngleHistogram[8], uint64_t(6));
        }

        void sharedFaceAndQuadraticNodes()
        {
            auto input = regular();
            auto output = input;
            const auto top = output.getNodes()[3];
            output.addNode(top.x, top.y, -top.z);
            output.addTetrahedron(0, 2, 1, 4, 1);
            // Face counts come from connectivity, independently of neighbor lists.
            const auto linear = RTetraMeshStatistics::measure(input, output);
            QCOMPARE(linear.faces, size_t(7));
            QCOMPARE(linear.exteriorFaces, size_t(6));
            QCOMPARE(linear.validElements, size_t(2));
            QCOMPARE(total(linear.faceAngleHistogram), uint64_t(24));
            QCOMPARE(total(linear.dihedralAngleHistogram), uint64_t(12));
            QVERIFY(output.createMidEdgeNodes() > 0);
            const auto quadratic = RTetraMeshStatistics::measure(input, output);
            QCOMPARE(quadratic.faces, linear.faces);
            QCOMPARE(quadratic.exteriorFaces, linear.exteriorFaces);
            QCOMPARE(quadratic.volume.minimum, linear.volume.minimum);
            QCOMPARE(quadratic.volume.maximum, linear.volume.maximum);
            QCOMPARE(quadratic.aspectHistogram, linear.aspectHistogram);
            QCOMPARE(quadratic.faceAngleHistogram, linear.faceAngleHistogram);
            QCOMPARE(quadratic.dihedralAngleHistogram, linear.dihedralAngleHistogram);
        }

        void emptyAndInvalid()
        {
            const RTetraMeshData empty;
            const auto none = RTetraMeshStatistics::measure(empty, empty);
            QCOMPARE(none.validElements, size_t(0));
            QCOMPARE(none.faces, size_t(0));
            QCOMPARE(total(none.aspectHistogram), uint64_t(0));
            QVERIFY(!std::isfinite(none.volume.minimum));
            auto data = regular();
            data.addTetrahedron(0, 2, 1, 3);  // Inverted.
            data.addTetrahedron(0, 1, 2, 2);  // Flat.
            data.addTetrahedron(0, 1, 2, 99); // Missing node.
            const auto result = RTetraMeshStatistics::measure(empty, data);
            QCOMPARE(result.validElements, size_t(1));
            QCOMPARE(result.invalidElements, size_t(3));
            QCOMPARE(total(result.aspectHistogram), uint64_t(1));
            QCOMPARE(total(result.faceAngleHistogram), uint64_t(12));
            QCOMPARE(total(result.dihedralAngleHistogram), uint64_t(6));
        }

        void nearlyFlatTetrahedron()
        {
            RTetraMeshData data;
            data.addNode(0, 0, 0);
            data.addNode(1, 0, 0);
            data.addNode(0, 1, 0);
            data.addNode(.25, .25, 1e-12);
            data.addTetrahedron(0, 1, 2, 3);
            const auto result = RTetraMeshStatistics::measure(RTetraMeshData(), data);
            QCOMPARE(result.validElements, size_t(1));
            QCOMPARE(result.invalidElements, size_t(0));
            QVERIFY(std::fabs(result.aspectRatio.minimum * 1e-12 - std::sqrt(2.0)) < 1e-12);
            QCOMPARE(result.aspectHistogram.back(), uint64_t(1));
            QVERIFY(result.dihedralAngle.minimum > 0.0);
            QVERIFY(result.dihedralAngle.minimum < 1e-8);
            QVERIFY(result.dihedralAngle.maximum > 179.9999);
            QCOMPARE(total(result.faceAngleHistogram), uint64_t(12));
            QCOMPARE(total(result.dihedralAngleHistogram), uint64_t(6));
        }

        void inputWinding()
        {
            auto input = regular();
            const auto output = input;
            std::swap(input.getTriangles()[0].nodeIDs[0], input.getTriangles()[0].nodeIDs[1]);
            std::swap(input.getTriangles()[2].nodeIDs[0], input.getTriangles()[2].nodeIDs[1]);
            const auto result = RTetraMeshStatistics::measure(input, output);
            QVERIFY(std::fabs(result.minimumInputFaceAngle - 60) < 1e-10);
            QVERIFY(std::fabs(result.minimumInputDihedralAngle - 70.5287793655093) < 1e-10);
        }

        void alignedLogRows()
        {
            static QStringList messages;
            messages.clear();
            auto &logger = RLogger::getInstance();
            const auto handler = logger.getLogHandler();
            const auto level = logger.getLevel();
            const auto restore = qScopeGuard([&]()
                                             {
                                                 logger.setLogHandler(handler);
                                                 logger.setLevel(level);
                                             });
            logger.setLevel(R_LOG_LEVEL_DETAIL);
            logger.setLogHandler([](const RMessage &message)
                                 {
                                     messages.append(message);
                                 });
            const auto data = regular();
            RTetraMeshStatistics::Counts counts;
            counts.inputSegments = counts.meshSegments = 6;
            counts.meshRegions = 1;
            const auto statistics = RTetraMeshStatistics::measure(data, data);
            const uint initialIndent = logger.getIndentLevel();
            RTetraMeshStatistics::print(data, data, counts, statistics);
            QCOMPARE(logger.getIndentLevel(), initialIndent);
            const QString log = messages.join("");
            QVERIFY(log.contains("Input mesh statistics"));
            QVERIFY(log.contains("Generated mesh statistics"));
            QVERIFY(log.contains("Mesh quality statistics"));
            QVERIFY(log.contains("Aspect ratio histogram"));
            QVERIFY(log.contains("Face angle histogram"));
            QVERIFY(log.contains("Dihedral angle histogram"));
            QVERIFY(log.contains("Minimum input dihedral angle ="));
            QVERIFY(log.contains("Face angle [deg]"));
            QVERIFY(!log.contains("asp.ratio"));
            QVERIFY(!log.contains("facangle"));
            QVERIFY(!log.contains('%'));
            QVERIFY(QRegularExpression("\\|\\s*Nodes\\s*\\|\\s*4\\s*\\|").match(log).hasMatch());
            int tableWidth = -1, rows = 0;
            for (const QString &line : log.split('\n'))
            {
                const QString trimmed = line.trimmed();
                if (trimmed.startsWith('+'))
                {
                    tableWidth = trimmed.size();
                }
                else if (trimmed.startsWith('|'))
                {
                    QVERIFY(tableWidth > 0);
                    QCOMPARE(trimmed.size(), tableWidth);
                    rows++;
                }
            }
            QVERIFY(rows > 60);
            for (const QString &message : messages)
            {
                if (message.contains('|') || message.contains("Aspect ratio ="))
                {
                    QCOMPARE(message.count('\n'), 1);
                }
            }
        }
};

QTEST_APPLESS_MAIN(TestTetraMeshStatistics)

#include "tst_tetra_mesh_statistics.moc"
