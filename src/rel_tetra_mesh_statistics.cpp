#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <QString>
#include <QScopeGuard>

#include <rbl_logger.h>

#include "rel_geometric_predicates.h"
#include "rel_tetra_mesh.h"
#include "rel_tetra_mesh_statistics.h"

namespace
{

    using Vector = std::array<double, 3>;
    using Face = std::array<int, 3>;
    using Edge = std::array<int, 2>;

    const std::array<double, 11> aspectLimits = {1.5, 2, 2.5, 3, 4, 6, 10, 15, 25, 50, 100};
    const std::array<double, 17> dihedralLimits = {5, 10, 20, 30, 40, 50, 60, 70, 80, 110, 120, 130, 140, 150, 160, 170, 175};

    Vector point(const RTetraMeshData &data, int id)
    {
        const auto &node = data.getNodes()[size_t(id)];
        return {node.x, node.y, node.z};
    }

    Vector subtract(const Vector &a, const Vector &b)
    {
        return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    }

    Vector cross(const Vector &a, const Vector &b)
    {
        return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
    }

    double dot(const Vector &a, const Vector &b)
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    double length(const Vector &v)
    {
        return std::hypot(v[0], v[1], v[2]);
    }

    bool normalize(Vector &v)
    {
        const double norm = length(v);
        if (!(norm > 0.0) || !std::isfinite(norm))
        {
            return false;
        }
        for (double &value : v)
        {
            value /= norm;
        }
        return true;
    }

    double angle(const Vector &a, const Vector &b)
    {
        static const double degrees = 180.0 / std::acos(-1.0);
        return std::atan2(length(cross(a, b)), dot(a, b)) * degrees;
    }

    template <size_t N>
    bool validNodes(const RTetraMeshData &data, const int (&ids)[N])
    {
        for (int id : ids)
        {
            if (id < 0 || size_t(id) >= data.getNodes().size())
            {
                return false;
            }
            const auto &node = data.getNodes()[size_t(id)];
            if (!std::isfinite(node.x) || !std::isfinite(node.y) || !std::isfinite(node.z))
            {
                return false;
            }
        }
        return true;
    }

    void measureInputAngles(const RTetraMeshData &input, RTetraMeshStatistics::Statistics &result)
    {
        std::map<Edge, std::vector<Vector>> normals;
        for (const auto &triangle : input.getTriangles())
        {
            if (!validNodes(input, triangle.nodeIDs))
            {
                continue;
            }
            for (int i = 0; i < 3; i++)
            {
                const int a = triangle.nodeIDs[i];
                const int b = triangle.nodeIDs[(i + 1) % 3];
                const int c = triangle.nodeIDs[(i + 2) % 3];
                Vector ab = subtract(point(input, b), point(input, a));
                Vector ac = subtract(point(input, c), point(input, a));
                if (normalize(ab) && normalize(ac))
                {
                    result.minimumInputFaceAngle = std::min(result.minimumInputFaceAngle, angle(ab, ac));
                }
                Edge key = {std::min(a, b), std::max(a, b)};
                // Orient every edge the same way; the angle between these normals
                // is the angle between the incident facets, regardless of winding.
                Vector edge = subtract(point(input, key[1]), point(input, key[0]));
                Vector toOpposite = subtract(point(input, c), point(input, key[0]));
                if (!normalize(edge) || !normalize(toOpposite))
                {
                    continue;
                }
                Vector normal = cross(edge, toOpposite);
                if (normalize(normal))
                {
                    normals[key].push_back(normal);
                }
            }
        }
        for (const auto &entry : normals)
        {
            const auto &faces = entry.second;
            for (size_t i = 0; i < faces.size(); i++)
            {
                for (size_t j = i + 1; j < faces.size(); j++)
                {
                    result.minimumInputDihedralAngle = std::min(result.minimumInputDihedralAngle, angle(faces[i], faces[j]));
                }
            }
        }
    }

    std::string number(double value)
    {
        return std::isfinite(value) ? QString::number(value, 'g', 7).toStdString() : "n/a";
    }

    using TableRow = std::vector<QString>;

    void printTable(const TableRow &header, const std::vector<TableRow> &rows)
    {
        std::vector<int> widths(header.size());
        for (size_t column = 0; column < header.size(); column++)
        {
            widths[column] = std::max(25, int(header[column].size()));
            for (const auto &row : rows)
            {
                widths[column] = std::max(widths[column], int(row[column].size()));
            }
        }
        QString border = "+";
        for (int width : widths)
        {
            border += " " + QString(width, '-') + " +";
        }
        auto printRow = [&](const TableRow &row)
        {
            QString line = "|";
            for (size_t column = 0; column < row.size(); column++)
            {
                line += QString(" %1 |").arg(row[column], column == 0 ? -widths[column] : widths[column]);
            }
            RLogger::info("%s\n", line.toUtf8().constData());
        };
        RLogger::info("%s\n", border.toUtf8().constData());
        printRow(header);
        RLogger::info("%s\n", border.toUtf8().constData());
        for (const auto &row : rows)
        {
            printRow(row);
        }
        RLogger::info("%s\n", border.toUtf8().constData());
    }

    template <typename Content>
    void printSection(const char *title, Content content)
    {
        RLogger::info("%s\n", title);
        RLogger::indent();
        const auto finish = qScopeGuard([]()
                                        {
                                            RLogger::unindent(false);
                                        });
        content();
    }

    template <size_t N>
    void printHistogram(const char *title, const char *quantity, const char *sample,
                        const std::array<uint64_t, N> &histogram, const std::array<QString, N> &labels)
    {
        printSection(title, [&]()
                     {
                         std::vector<TableRow> rows;
                         for (size_t i = 0; i < N; i++)
                         {
                             rows.push_back({labels[i], QString::number(qulonglong(histogram[i]))});
                         }
                         printTable({quantity, sample}, rows);
                     });
    }

}

void RTetraMeshStatistics::Range::add(double value)
{
    if (std::isfinite(value))
    {
        minimum = std::min(minimum, value);
        maximum = std::max(maximum, value);
    }
}

RTetraMeshStatistics::Statistics RTetraMeshStatistics::measure(const RTetraMeshData &input, const RTetraMeshData &output)
{
    Statistics result;
    std::map<Face, size_t> faces;
    for (const auto &tetrahedron : output.getTetrahedra())
    {
        if (!validNodes(output, tetrahedron.nodeIDs))
        {
            result.invalidElements++;
            continue;
        }
        Vector points[4];
        for (int i = 0; i < 4; i++)
        {
            points[i] = point(output, tetrahedron.nodeIDs[i]);
            const int *f = RTetraMesh::faceVertices[i];
            Face key = {tetrahedron.nodeIDs[f[0]], tetrahedron.nodeIDs[f[1]], tetrahedron.nodeIDs[f[2]]};
            std::sort(key.begin(), key.end());
            faces[key]++;
        }
        double longest = 0.0;
        for (int i = 0; i < 4; i++)
        {
            for (int j = i + 1; j < 4; j++)
            {
                longest = std::max(longest, length(subtract(points[i], points[j])));
            }
        }
        if (!(longest > 0.0) || !std::isfinite(longest))
        {
            result.invalidElements++;
            continue;
        }
        Vector scaled[4] = {};
        for (int i = 1; i < 4; i++)
        {
            scaled[i] = subtract(points[i], points[0]);
            for (double &value : scaled[i])
            {
                value /= longest;
            }
        }
        const double determinant = RGeometricPredicates::orient3d(scaled[0].data(), scaled[1].data(), scaled[2].data(), scaled[3].data());
        const double volume = determinant * longest * longest * longest / 6.0;
        if (!(determinant > 0.0) || !std::isfinite(volume) || !(volume > 0.0))
        {
            result.invalidElements++;
            continue;
        }
        Vector normals[4];
        double largestNormal = 0.0;
        bool valid = true;
        for (int i = 0; i < 4; i++)
        {
            const int *f = RTetraMesh::faceVertices[i];
            normals[i] = cross(subtract(scaled[f[1]], scaled[f[0]]), subtract(scaled[f[2]], scaled[f[0]]));
            largestNormal = std::max(largestNormal, length(normals[i]));
            valid = normalize(normals[i]) && valid;
        }
        const double aspect = largestNormal / determinant;
        if (!valid || !std::isfinite(aspect))
        {
            result.invalidElements++;
            continue;
        }
        result.validElements++;
        result.volume.add(volume);
        result.aspectRatio.add(aspect);
        // TetGen puts an exact upper-bound value in the lower aspect-ratio bin.
        const size_t aspectBin = size_t(std::lower_bound(aspectLimits.begin(), aspectLimits.end(), aspect) - aspectLimits.begin());
        result.aspectHistogram[aspectBin]++;
        for (int i = 0; i < 4; i++)
        {
            for (int j = i + 1; j < 4; j++)
            {
                result.edge.add(length(subtract(points[i], points[j])));
                Vector opposite = normals[j];
                for (double &value : opposite)
                {
                    value = -value;
                }
                const double dihedral = angle(normals[i], opposite);
                result.dihedralAngle.add(dihedral);
                const size_t bin = size_t(std::upper_bound(dihedralLimits.begin(), dihedralLimits.end(), dihedral) - dihedralLimits.begin());
                result.dihedralAngleHistogram[bin]++;
            }
            const int *f = RTetraMesh::faceVertices[i];
            for (int j = 0; j < 3; j++)
            {
                Vector a = subtract(scaled[f[(j + 1) % 3]], scaled[f[j]]);
                Vector b = subtract(scaled[f[(j + 2) % 3]], scaled[f[j]]);
                normalize(a);
                normalize(b);
                const double faceAngle = angle(a, b);
                result.faceAngle.add(faceAngle);
                const size_t bin = std::min(size_t(faceAngle / 10.0), size_t(17));
                result.faceAngleHistogram[bin]++;
            }
        }
    }
    result.faces = faces.size();
    for (const auto &face : faces)
    {
        result.exteriorFaces += face.second == 1;
    }
    measureInputAngles(input, result);
    return result;
}

void RTetraMeshStatistics::print(const RTetraMeshData &input, const RTetraMeshData &output,
                                 const Counts &counts, const Statistics &statistics)
{
    auto count = [](uint64_t value)
    {
        return QString::number(qulonglong(value));
    };
    printSection("Input mesh statistics", [&]()
                 {
                     printTable({"Entity", "Count"}, {{"Nodes", count(input.getNodes().size())},
                                                      {"Surface triangles", count(input.getTriangles().size())},
                                                      {"Constrained edges", count(counts.inputSegments)},
                                                      // RTetraMeshData does not contain explicit hole seed points.
                                                      {"Hole seeds", "0"},
                                                      {"Regions", count(counts.inputRegions)}});
                 });
    printSection("Generated mesh statistics", [&]()
                 {
                     std::vector<TableRow> rows = {
                         {"Nodes", count(output.getNodes().size())},
                         {"Tetrahedra", count(output.getTetrahedra().size())},
                         {"Faces", count(statistics.faces)},
                         {"Exterior boundary faces", count(statistics.exteriorFaces)},
                         {"Faces on input facets", count(output.getTriangles().size())},
                         {"Edges on input segments", count(counts.meshSegments)},
                         {"Regions", count(counts.meshRegions)},
                         {"Added nodes on facets", count(counts.steinerOnFacets)},
                         {"Added nodes on segments", count(counts.steinerOnSegments)},
                         {"Added interior nodes", count(counts.steinerInside)}};
                     if (counts.midEdgeNodes > 0)
                     {
                         rows.push_back({"Mid-edge nodes", count(counts.midEdgeNodes)});
                     }
                     if (statistics.invalidElements > 0)
                     {
                         rows.push_back({"Invalid tetrahedra", count(statistics.invalidElements)});
                     }
                     printTable({"Entity", "Count"}, rows);
                 });
    printSection("Mesh quality statistics", [&]()
                 {
                     auto rangeRow = [](const char *name, const Range &range)
                     {
                         return TableRow{name, QString::fromStdString(number(range.minimum)), QString::fromStdString(number(range.maximum))};
                     };
                     printTable({"Quantity", "Minimum", "Maximum"}, {rangeRow("Volume", statistics.volume),
                                                                     rangeRow("Edge length", statistics.edge),
                                                                     rangeRow("Aspect ratio", statistics.aspectRatio),
                                                                     rangeRow("Face angle [deg]", statistics.faceAngle),
                                                                     rangeRow("Dihedral angle [deg]", statistics.dihedralAngle)});
                     RLogger::info("Aspect ratio = longest edge length / minimum tetrahedron altitude.\n");
                     if (std::isfinite(statistics.minimumInputFaceAngle))
                     {
                         RLogger::info("Minimum input face angle = %.7g [deg]\n", statistics.minimumInputFaceAngle);
                     }
                     if (std::isfinite(statistics.minimumInputDihedralAngle))
                     {
                         RLogger::info("Minimum input dihedral angle = %.7g [deg]\n", statistics.minimumInputDihedralAngle);
                     }
                     if (statistics.invalidElements > 0)
                     {
                         RLogger::info("Quality ranges and histograms exclude invalid tetrahedra.\n");
                     }
                 });

    std::array<QString, 12> aspectLabels;
    for (size_t i = 0; i < aspectLabels.size(); i++)
    {
        if (i == 0)
        {
            aspectLabels[i] = "<= 1.5";
        }
        else if (i == aspectLimits.size())
        {
            aspectLabels[i] = "> 100";
        }
        else
        {
            aspectLabels[i] = QString("%1 - %2").arg(aspectLimits[i - 1]).arg(aspectLimits[i]);
        }
    }
    printHistogram("Aspect ratio histogram", "Aspect ratio", "Elements", statistics.aspectHistogram, aspectLabels);
    std::array<QString, 18> faceLabels, dihedralLabels;
    for (size_t i = 0; i < faceLabels.size(); i++)
    {
        faceLabels[i] = QString("%1 - %2").arg(i * 10).arg((i + 1) * 10);
        const double lower = i == 0 ? 0.0 : dihedralLimits[i - 1];
        const double upper = i == dihedralLimits.size() ? 180.0 : dihedralLimits[i];
        dihedralLabels[i] = QString("%1 - %2").arg(lower).arg(upper);
    }
    printHistogram("Face angle histogram", "Face angle [deg]", "Angles", statistics.faceAngleHistogram, faceLabels);
    printHistogram("Dihedral angle histogram", "Dihedral angle [deg]", "Angles", statistics.dihedralAngleHistogram, dihedralLabels);
}
