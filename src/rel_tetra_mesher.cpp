#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <vector>

#include <rbl_error.h>
#include <rbl_logger.h>
#include <rbl_progress.h>

#include "rel_boundary_recovery.h"
#include "rel_delaunay_tetrahedralizer.h"
#include "rel_geometric_predicates.h"
#include "rel_tetra_mesh.h"
#include "rel_tetra_mesh_refiner.h"
#include "rel_tetra_mesher.h"
#include "rel_tetra_mesh_optimizer.h"
#include "rel_tetra_mesh_quality.h"
#include "rel_tetra_mesh_statistics.h"

namespace
{

    //! Length below which nothing is split, relative to size of the model.
    const double minLengthRatio = 1.0e-8;
    //! Ratio of length of a boundary edge to size of the model, below which added vertices at it are removed after recovery.
    const double cleanUpLengthRatio = 1.0e-4;
    //! Number of nodes which can be added if it is not limited by options.
    const uint64_t defaultMaxSteinerNodes = 5000000;
    //! Number of tetrahedra of each region which are used to find a region of input mesh.
    const int nRegionSamples = 5;

    void validateInput(const RTetraMeshData &input, bool useTetrahedra)
    {
        const int nNodes = int(input.getNodes().size());

        for (const RTetraMeshData::Node &rNode : input.getNodes())
        {
            if (!std::isfinite(rNode.x) || !std::isfinite(rNode.y) || !std::isfinite(rNode.z))
            {
                throw RError(RError::Type::InvalidInput, R_ERROR_REF, "Node coordinates are not finite numbers.");
            }
        }
        for (const RTetraMeshData::Edge &rEdge : input.getEdges())
        {
            for (int nodeID : rEdge.nodeIDs)
            {
                if (nodeID < 0 || nodeID >= nNodes)
                {
                    throw RError(RError::Type::InvalidInput, R_ERROR_REF, "Edge refers to a node which does not exist (%d).", nodeID);
                }
            }
        }
        for (const RTetraMeshData::Triangle &rTriangle : input.getTriangles())
        {
            for (int nodeID : rTriangle.nodeIDs)
            {
                if (nodeID < 0 || nodeID >= nNodes)
                {
                    throw RError(RError::Type::InvalidInput, R_ERROR_REF, "Triangle refers to a node which does not exist (%d).", nodeID);
                }
            }
        }
        if (useTetrahedra)
        {
            for (const RTetraMeshData::Tetrahedron &rTetrahedron : input.getTetrahedra())
            {
                for (int nodeID : rTetrahedron.nodeIDs)
                {
                    if (nodeID < 0 || nodeID >= nNodes)
                    {
                        throw RError(RError::Type::InvalidInput, R_ERROR_REF, "Tetrahedron refers to a node which does not exist (%d).", nodeID);
                    }
                }
            }
        }
    }

    //! Find faces of input tetrahedra which lie on a boundary of a region.
    std::vector<std::array<int, 3>> findRegionFaces(const RTetraMeshData &input, const std::vector<int> &nodeMap)
    {
        static const int faceNodes[4][3] = {{1, 3, 2}, {0, 2, 3}, {0, 3, 1}, {0, 1, 2}};

        struct FaceInfo
        {
            std::array<int, 3> nodeIDs;
            int marker;
            int count;
            bool boundary;
        };

        std::map<std::array<int, 3>, FaceInfo> faces;

        for (const RTetraMeshData::Tetrahedron &rTetrahedron : input.getTetrahedra())
        {
            for (int i = 0; i < 4; i++)
            {
                const std::array<int, 3> nodeIDs = {nodeMap[size_t(rTetrahedron.nodeIDs[faceNodes[i][0]])],
                                                    nodeMap[size_t(rTetrahedron.nodeIDs[faceNodes[i][1]])],
                                                    nodeMap[size_t(rTetrahedron.nodeIDs[faceNodes[i][2]])]};
                std::array<int, 3> key = nodeIDs;
                std::sort(key.begin(), key.end());

                auto iter = faces.find(key);
                if (iter == faces.end())
                {
                    faces[key] = FaceInfo{nodeIDs, rTetrahedron.marker, 1, false};
                }
                else
                {
                    iter->second.count++;
                    if (iter->second.marker != rTetrahedron.marker)
                    {
                        iter->second.boundary = true;
                    }
                }
            }
        }

        std::vector<std::array<int, 3>> regionFaces;
        for (const auto &rItem : faces)
        {
            if (rItem.second.count == 1 || rItem.second.boundary)
            {
                regionFaces.push_back(rItem.second.nodeIDs);
            }
        }
        return regionFaces;
    }

    //! Find marker of an input tetrahedron which contains a point.
    bool findInputMarker(const RTetraMeshData &input, const double *p, int &marker)
    {
        const std::vector<RTetraMeshData::Node> &rNodes = input.getNodes();
        for (const RTetraMeshData::Tetrahedron &rTetrahedron : input.getTetrahedra())
        {
            double x[4][3];
            bool outside = false;
            for (int j = 0; j < 3 && !outside; j++)
            {
                double lower = std::numeric_limits<double>::max();
                double upper = std::numeric_limits<double>::lowest();
                for (int i = 0; i < 4; i++)
                {
                    const RTetraMeshData::Node &rNode = rNodes[size_t(rTetrahedron.nodeIDs[i])];
                    x[i][j] = j == 0 ? rNode.x : (j == 1 ? rNode.y : rNode.z);
                    lower = std::min(lower, x[i][j]);
                    upper = std::max(upper, x[i][j]);
                }
                outside = p[j] < lower || p[j] > upper;
            }
            if (outside)
            {
                continue;
            }
            const double orientation = RGeometricPredicates::orient3d(x[0], x[1], x[2], x[3]);
            if (orientation == 0.0)
            {
                continue;
            }
            const double sign = orientation > 0.0 ? 1.0 : -1.0;
            if (sign * RGeometricPredicates::orient3d(x[1], x[3], x[2], p) >= 0.0 && sign * RGeometricPredicates::orient3d(x[0], x[2], x[3], p) >= 0.0 && sign * RGeometricPredicates::orient3d(x[0], x[3], x[1], p) >= 0.0 && sign * RGeometricPredicates::orient3d(x[0], x[1], x[2], p) >= 0.0)
            {
                marker = rTetrahedron.marker;
                return true;
            }
        }
        return false;
    }

}

void RTetraMesher::generate(const RTetraMeshData &input, const RTetraMeshOptions &options, RTetraMeshData &output)
{
    generate(input, options, output, nullptr);
}

void RTetraMesher::generate(const RTetraMeshData &input, const RTetraMeshOptions &options, RTetraMeshData &output, RTetraMeshProfile *resultProfile)
{
    RLogger::info("Preparing volume mesh\n");
    RLogger::indent();

    RTetraMeshProfile profile;
    RTetraMeshProfile::Timer totalTimer(profile.total);
    const std::vector<RTetraMeshData::Node> &rInputNodes = input.getNodes();
    const int nInputNodes = int(rInputNodes.size());
    const bool reconstruct = options.getReconstruct() && !input.getTetrahedra().empty();

    validateInput(input, reconstruct);

    if (options.getElementOrder() != 1 && options.getElementOrder() != 2)
    {
        RLogger::unindent();
        throw RError(RError::Type::InvalidInput, R_ERROR_REF, "Order of generated elements must be 1 or 2 (%d).", options.getElementOrder());
    }

    if (input.getTriangles().empty() && !reconstruct)
    {
        RLogger::unindent();
        throw RError(RError::Type::InvalidInput, R_ERROR_REF, "There are no triangles to generate tetrahedra from.");
    }

    // Delaunay tetrahedralization of input nodes.

    RTetraMesh mesh;
    mesh.setCoplanarTolerance(std::max(options.getTolerance(), 1.0e-12));
    double lower[3] = {0.0, 0.0, 0.0};
    double upper[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < nInputNodes; i++)
    {
        const double x[3] = {rInputNodes[size_t(i)].x, rInputNodes[size_t(i)].y, rInputNodes[size_t(i)].z};
        mesh.addPoint(x[0], x[1], x[2]);
        for (int j = 0; j < 3; j++)
        {
            lower[j] = i == 0 ? x[j] : std::min(lower[j], x[j]);
            upper[j] = i == 0 ? x[j] : std::max(upper[j], x[j]);
        }
    }
    const double modelSize = std::hypot(upper[0] - lower[0], upper[1] - lower[1], upper[2] - lower[2]);

    RLogger::info("Tetrahedralizing nodes\n");
    if (options.getVerbose())
    {
        RLogger::info("Tetrahedralizing %d nodes.\n", nInputNodes);
    }

    std::vector<int> duplicates;
    {
        RTetraMeshProfile::Timer timer(profile.tetrahedralization);
        RDelaunayTetrahedralizer::tetrahedralize(mesh, duplicates);
    }

    std::vector<int> nodeMap(static_cast<size_t>(nInputNodes));
    int nDuplicates = 0;
    for (int i = 0; i < nInputNodes; i++)
    {
        nodeMap[size_t(i)] = duplicates[size_t(i)] >= 0 ? duplicates[size_t(i)] : i;
        if (duplicates[size_t(i)] >= 0)
        {
            nDuplicates++;
        }
    }
    if (nDuplicates > 0)
    {
        RLogger::warning("Found %d nodes with the same or nearly the same coordinates as other nodes, they are not used in generated mesh.\n", nDuplicates);
    }

    // Boundary.

    RLogger::info("Recovering mesh boundary\n");
    RTetraMeshProfile::Timer boundaryTimer(profile.boundary);
    RBoundaryRecovery boundary(mesh, options.getTolerance());
    boundary.setMinLength(minLengthRatio * modelSize);

    if (options.getUseSizeFunction())
    {
        for (int i = 0; i < nInputNodes; i++)
        {
            const double size = rInputNodes[size_t(i)].size;
            if (size > 0.0 && boundary.getVertexSize(nodeMap[size_t(i)]) <= 0.0)
            {
                boundary.setVertexSize(nodeMap[size_t(i)], size);
            }
        }
    }

    for (const RTetraMeshData::Triangle &rTriangle : input.getTriangles())
    {
        boundary.addTriangle(nodeMap[size_t(rTriangle.nodeIDs[0])],
                             nodeMap[size_t(rTriangle.nodeIDs[1])],
                             nodeMap[size_t(rTriangle.nodeIDs[2])],
                             rTriangle.marker);
    }
    if (reconstruct)
    {
        // Regions of input mesh do not have to be separated by triangles.
        for (const std::array<int, 3> &rFace : findRegionFaces(input, nodeMap))
        {
            boundary.addTriangle(rFace[0], rFace[1], rFace[2], RBoundaryRecovery::noMarker);
        }
    }
    for (const RTetraMeshData::Edge &rEdge : input.getEdges())
    {
        boundary.addEdge(nodeMap[size_t(rEdge.nodeIDs[0])], nodeMap[size_t(rEdge.nodeIDs[1])], rEdge.marker);
    }

    if (options.getVerbose())
    {
        RLogger::info("Recovering boundary.\n");
    }

    const uint64_t maxSteinerNodes = options.getMaxSteinerNodes() > 0 ? uint64_t(options.getMaxSteinerNodes()) : defaultMaxSteinerNodes;

    // Recovery which does not converge is stopped by a limit of its own, unless the limit given by options is lower.
    boundary.setMaxSteinerVertices(std::min(maxSteinerNodes, 100 * uint64_t(nInputNodes) + 100000));
    boundary.initialize();
    const size_t nInputSegments = size_t(std::count_if(boundary.getSubSegments().begin(), boundary.getSubSegments().end(),
                                                       [](const RBoundaryRecovery::SubSegment &segment)
                                                       {
                                                           return segment.alive;
                                                       }));
    {
        // Recovery by splitting fails when facets are triangulated differently in the tetrahedralization:
        // parts of their missing subfaces converge to a point. Only then is recovery started again,
        // allowing facets to take edges of the tetrahedralization, which leaves a worse boundary.
        RTetraMesh initialMesh = mesh;
        RBoundaryRecovery initialBoundary(initialMesh, boundary);
        try
        {
            boundary.recover();
        }
        catch (const RError &rError)
        {
            RLogger::warning("Recovering mesh boundary again, with facets taking edges of the tetrahedralization: %s\n",
                             rError.getMessage().toUtf8().constData());
            mesh = std::move(initialMesh);
            boundary.restore(initialBoundary);
            boundary.setFacetEdgeInsertion(true);
            boundary.recover();
        }
    }
    {
        // Facets which are triangulated differently in the tetrahedralization leave clusters of
        // added vertices and needle-like subfaces, which refinement can not handle. Those vertices
        // are removed and facets are made Delaunay again, as far as the tetrahedralization allows.
        const uint64_t nRemoved = boundary.removeShortEdges(cleanUpLengthRatio * modelSize);
        boundary.recover();
        const uint64_t nFlips = boundary.restoreFacetDelaunay();
        if (options.getVerbose())
        {
            RLogger::info("Removed %llu added boundary nodes and flipped %llu facet edges.\n",
                          static_cast<unsigned long long>(nRemoved), static_cast<unsigned long long>(nFlips));
        }
    }
    boundary.setFacetEdgeInsertion(false);

    if (options.getVerbose())
    {
        RLogger::info("Boundary was recovered with %llu nodes added on edges and %llu nodes added on triangles.\n",
                      static_cast<unsigned long long>(boundary.getNSteinerVertices(RBoundaryRecovery::VertexSegment)),
                      static_cast<unsigned long long>(boundary.getNSteinerVertices(RBoundaryRecovery::VertexFacet)));
    }

    RLogger::info("Classifying mesh regions");
    int nRegions = boundary.classifyRegions();
    if (nRegions == 0)
    {
        throw RError(RError::Type::InvalidInput, R_ERROR_REF, "Triangles do not enclose any volume.");
    }

    boundaryTimer.finish();
    RTetraMeshOptimizer optimizer(mesh, boundary, options);
    {
        RTetraMeshProfile::Timer timer(profile.initialRepair);
        optimizer.optimize(false);
    }

    // Refinement.

    if (options.getQualityMesh())
    {
        RTetraMeshProfile::Timer refinementTimer(profile.refinement);

        // The limit is checked before a tetrahedron is refined, which may need more than one node.
        boundary.setMaxSteinerVertices(std::numeric_limits<uint64_t>::max());

        RTetraMeshRefiner refiner(mesh, boundary, options);
        refiner.setMinLength(minLengthRatio * modelSize);
        refiner.setMaxSteinerVertices(maxSteinerNodes);
        refiner.refine();
        profile.addWork(refiner.getProfile());
        if (options.getVerbose())
        {
            RLogger::info("Refinement safeguards rejected %llu interior insertions and %llu boundary splits; constrained cells may remain below requested quality.\n",
                          static_cast<unsigned long long>(refiner.getSpacingRejections()), static_cast<unsigned long long>(refiner.getBoundaryRejections()));
        }

        if (refiner.getStopped())
        {
            RLogger::warning("Mesh refinement was stopped after %llu nodes were added, some elements do not meet required quality.\n",
                             static_cast<unsigned long long>(boundary.getNSteinerVertices()));
        }

        RLogger::info("Recovering refined mesh boundary\n");
        boundary.recover();
        RLogger::info("Classifying refined mesh regions\n");
        nRegions = boundary.classifyRegions();
        if (nRegions == 0)
        {
            throw RError(RError::Type::Application, R_ERROR_REF, "Boundary was lost during mesh refinement.");
        }
    }

    if (options.getQualityMesh())
    {
        if (options.getVerbose())
        {
            RLogger::info("Optimizing mesh with bounded patch repairs.\n");
        }
        RTetraMeshProfile::Timer timer(profile.optimization);
        optimizer.optimize();
    }

    // Output.

    RLogger::info("Preparing generated mesh output\n");
    output.clear();

    std::vector<RTetraMeshData::Node> &rNodes = output.getNodes();
    rNodes = rInputNodes;
    // Added nodes which were removed again are left out.
    std::vector<int> pointMap(size_t(mesh.getNPoints()));
    for (int i = 0; i < mesh.getNPoints(); i++)
    {
        pointMap[size_t(i)] = i;
    }
    for (int i = nInputNodes; i < mesh.getNPoints(); i++)
    {
        if (!mesh.isPointUsed(i))
        {
            pointMap[size_t(i)] = -1;
            continue;
        }
        pointMap[size_t(i)] = int(rNodes.size());
        const double *p = mesh.getPoint(i);
        output.addNode(p[0], p[1], p[2], RTetraMeshData::noMarker, boundary.getVertexSize(i));
    }

    const int nSlots = mesh.getNTetraSlots();
    std::vector<int> tetraMap(static_cast<size_t>(nSlots), RTetraMeshData::noTetrahedron);
    int nTetrahedra = 0;
    for (int i = 0; i < nSlots; i++)
    {
        if (mesh.isAlive(i) && !mesh.isGhost(i) && mesh.getTetra(i).region >= 0)
        {
            tetraMap[size_t(i)] = nTetrahedra++;
        }
    }

    // Markers of regions.
    std::vector<int> regionMarkers(static_cast<size_t>(nRegions));
    for (int i = 0; i < nRegions; i++)
    {
        regionMarkers[size_t(i)] = i;
    }
    if (reconstruct)
    {
        int maxMarker = -1;
        for (const RTetraMeshData::Tetrahedron &rTetrahedron : input.getTetrahedra())
        {
            maxMarker = std::max(maxMarker, rTetrahedron.marker);
        }

        std::vector<std::map<int, int>> votes(static_cast<size_t>(nRegions));
        std::vector<int> nSamples(size_t(nRegions), 0);
        for (int i = 0; i < nSlots; i++)
        {
            if (tetraMap[size_t(i)] == RTetraMeshData::noTetrahedron)
            {
                continue;
            }
            const RTetraMesh::Tetra &rTetra = mesh.getTetra(i);
            if (nSamples[size_t(rTetra.region)] >= nRegionSamples)
            {
                continue;
            }
            nSamples[size_t(rTetra.region)]++;
            double center[3] = {0.0, 0.0, 0.0};
            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    center[k] += 0.25 * mesh.getPoint(rTetra.v[j])[k];
                }
            }
            int marker;
            if (findInputMarker(input, center, marker))
            {
                votes[size_t(rTetra.region)][marker]++;
            }
        }
        for (int i = 0; i < nRegions; i++)
        {
            int bestCount = 0;
            int bestMarker = RTetraMeshData::noMarker;
            for (const auto &rVote : votes[size_t(i)])
            {
                if (rVote.second > bestCount)
                {
                    bestCount = rVote.second;
                    bestMarker = rVote.first;
                }
            }
            // Region which was not a part of input mesh.
            regionMarkers[size_t(i)] = bestCount > 0 ? bestMarker : ++maxMarker;
        }
    }

    std::vector<RTetraMeshData::Tetrahedron> &rTetrahedra = output.getTetrahedra();
    rTetrahedra.reserve(size_t(nTetrahedra));
    for (int i = 0; i < nSlots; i++)
    {
        if (tetraMap[size_t(i)] == RTetraMeshData::noTetrahedron)
        {
            continue;
        }
        const RTetraMesh::Tetra &rTetra = mesh.getTetra(i);
        RTetraMeshData::Tetrahedron tetrahedron;
        for (int j = 0; j < 4; j++)
        {
            tetrahedron.nodeIDs[j] = pointMap[size_t(rTetra.v[j])];
            tetrahedron.neighborIDs[j] = tetraMap[size_t(rTetra.n[j])];
        }
        tetrahedron.marker = regionMarkers[size_t(rTetra.region)];
        rTetrahedra.push_back(tetrahedron);
    }

    for (const RBoundaryRecovery::SubFace &rFace : boundary.getSubFaces())
    {
        if (rFace.alive)
        {
            output.addTriangle(pointMap[size_t(rFace.v[0])], pointMap[size_t(rFace.v[1])], pointMap[size_t(rFace.v[2])], rFace.marker);
        }
    }
    for (const RBoundaryRecovery::SubSegment &rSegment : boundary.getSubSegments())
    {
        if (rSegment.alive && rSegment.marker != RBoundaryRecovery::noMarker)
        {
            output.addEdge(pointMap[size_t(rSegment.v[0])], pointMap[size_t(rSegment.v[1])], rSegment.marker);
        }
    }

    // Quadratic elements. Nodes are added after everything else was done on a linear mesh,
    // which they do not change: the nodes lie in the middle of straight edges.
    int nMidEdgeNodes = 0;
    if (options.getElementOrder() == 2)
    {
        RLogger::info("Creating mid-edge nodes\n");
        nMidEdgeNodes = output.createMidEdgeNodes();
    }

    profile.addWork(optimizer.getProfile());
    RLogger::info("Checking mesh quality\n");
    RTetraMeshQuality::Summary quality;
    RTetraMeshStatistics::Statistics statistics;
    {
        RTetraMeshProfile::Timer timer(profile.diagnostics);
        quality = RTetraMeshQuality::summarize(output, options);
        if (options.getVerbose())
        {
            statistics = RTetraMeshStatistics::measure(input, output);
        }
    }
    totalTimer.finish(); // Excludes optional report serialization and logging.
    if (!options.getQualityReportPath().empty())
    {
        RLogger::info("Writing mesh quality report\n");
        const RTetraMeshQuality::Work work = {optimizer.getQueueVisits(), optimizer.getBudgetStops(), optimizer.getInsertedPoints(), optimizer.getOrderedVisits(), optimizer.getBoundaryPoints(), &profile};
        RTetraMeshQuality::writeReport(output, options, optimizer.getPatchAttempts(), options.getQualityReportPath(), &work);
    }
    if (resultProfile)
    {
        *resultProfile = profile;
    }
    if (options.getVerbose())
    {
        RLogger::info("Reporting mesh statistics\n");
        RLogger::indent();
        RTetraMeshStatistics::Counts counts;
        counts.inputSegments = nInputSegments;
        counts.meshSegments = size_t(std::count_if(boundary.getSubSegments().begin(), boundary.getSubSegments().end(),
                                                   [](const RBoundaryRecovery::SubSegment &segment)
                                                   {
                                                       return segment.alive;
                                                   }));
        counts.meshRegions = size_t(nRegions);
        counts.steinerOnFacets = boundary.getNSteinerVertices(RBoundaryRecovery::VertexFacet);
        counts.steinerOnSegments = boundary.getNSteinerVertices(RBoundaryRecovery::VertexSegment);
        counts.steinerInside = boundary.getNSteinerVertices(RBoundaryRecovery::VertexVolume);
        counts.midEdgeNodes = size_t(nMidEdgeNodes);
        if (reconstruct)
        {
            std::set<int> regions;
            for (const auto &tetrahedron : input.getTetrahedra())
            {
                regions.insert(tetrahedron.marker);
            }
            counts.inputRegions = regions.size();
        }
        RTetraMeshStatistics::print(input, output, counts, statistics);

        int textWidth = 26;
        int numWidth = 25;
        RLogger::info("Mesh phases\n");
        RLogger::indent();
        RLogger::info("%-*s %*.3f [s]\n", textWidth, "tetrahedralization:", numWidth, profile.tetrahedralization);
        RLogger::info("%-*s %*.3f [s]\n", textWidth, "boundary:", numWidth, profile.boundary);
        RLogger::info("%-*s %*.3f [s]\n", textWidth, "initial repair:", numWidth, profile.initialRepair);
        RLogger::info("%-*s %*.3f [s]\n", textWidth, "refinement:", numWidth, profile.refinement);
        RLogger::info("%-*s %*.3f [s]\n", textWidth, "optimization:", numWidth, profile.optimization);
        RLogger::info("%-*s %*.3f [s]\n", textWidth, "diagnostics:", numWidth, profile.diagnostics);
        RLogger::info("%-*s %*.3f [s]\n", textWidth, "total:", numWidth, profile.total);
        RLogger::unindent(false);
        RLogger::info("Mesh work\n");
        RLogger::indent();
        RLogger::info("%-*s %*llu\n", textWidth, "quality measurements:", numWidth, static_cast<unsigned long long>(profile.qualityMeasurements));
        RLogger::info("%-*s %*llu\n", textWidth, "cache hits:", numWidth, static_cast<unsigned long long>(profile.qualityCacheHits));
        RLogger::info("%-*s %*llu\n", textWidth, "star searches:", numWidth, static_cast<unsigned long long>(profile.starQueries));
        RLogger::info("%-*s %*llu\n", textWidth, "star cache hits:", numWidth, static_cast<unsigned long long>(profile.starCacheHits));
        RLogger::info("%-*s %*llu\n", textWidth, "local trial replacements:", numWidth, static_cast<unsigned long long>(profile.trialReplacements));
        RLogger::unindent(false);
        RLogger::info("Volume optimization visited %llu local cell patches, accepted %llu off-centers and reached %llu queue budgets.\n",
                      static_cast<unsigned long long>(optimizer.getQueueVisits()), static_cast<unsigned long long>(optimizer.getInsertedPoints()), static_cast<unsigned long long>(optimizer.getBudgetStops()));
        RLogger::info("Sub-degree angles: %llu elements have internal face pairs, %llu have prescribed face pairs (counts may overlap).\n",
                      static_cast<unsigned long long>(quality.freeBelowOne), static_cast<unsigned long long>(quality.prescribedBelowOne));
        RLogger::info("Mesh quality\n");
        RLogger::indent();
        RLogger::info("%-*s %*.3g [deg]\n", textWidth, "minimum dihedral:", numWidth, quality.minAngle);
        RLogger::info("%-*s %*.3g\n", textWidth, "radius/edge p95:", numWidth, quality.radiusEdge95);
        RLogger::info("%-*s %*.3g\n", textWidth, "maximum:", numWidth, quality.maxRadiusEdge);
        RLogger::info("%-*s %*llu\n", textWidth, "nearly flat:", numWidth, static_cast<unsigned long long>(quality.nearFlat));
        RLogger::unindent(false);
        RLogger::unindent(false);
    }
    if (options.getQualityMesh() && !quality.targetsMet())
    {
        RLogger::warning("Mesh quality targets remain unmet: %llu radius/edge, %llu dihedral, %llu volume/size, %llu nearly flat, %llu invalid elements (%llu elements have geometry-constrained angles).\n",
                         static_cast<unsigned long long>(quality.badRadius), static_cast<unsigned long long>(quality.badAngle),
                         static_cast<unsigned long long>(quality.tooLarge), static_cast<unsigned long long>(quality.nearFlat),
                         static_cast<unsigned long long>(quality.invalid), static_cast<unsigned long long>(quality.constrainedAngle));
    }

    if (options.getVerbose())
    {
        RLogger::info("Generated %d tetrahedra in %d volumes with %llu nodes added (%llu on edges, %llu on triangles, %llu inside).\n",
                      nTetrahedra,
                      nRegions,
                      static_cast<unsigned long long>(boundary.getNSteinerVertices()),
                      static_cast<unsigned long long>(boundary.getNSteinerVertices(RBoundaryRecovery::VertexSegment)),
                      static_cast<unsigned long long>(boundary.getNSteinerVertices(RBoundaryRecovery::VertexFacet)),
                      static_cast<unsigned long long>(boundary.getNSteinerVertices(RBoundaryRecovery::VertexVolume)));
        if (options.getElementOrder() == 2)
        {
            RLogger::info("Elements are quadratic, %d nodes were added in the middle of edges.\n", nMidEdgeNodes);
        }
    }
    RLogger::unindent();
}
