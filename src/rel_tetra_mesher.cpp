#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <vector>

#include <rbl_error.h>
#include <rbl_logger.h>

#include "rel_boundary_recovery.h"
#include "rel_delaunay_tetrahedralizer.h"
#include "rel_geometric_predicates.h"
#include "rel_tetra_mesh.h"
#include "rel_tetra_mesh_progress.h"
#include "rel_tetra_mesh_refiner.h"
#include "rel_tetra_mesher.h"
#include "rel_tetra_mesh_optimizer.h"
#include "rel_tetra_mesh_quality.h"

namespace
{

    //! Length below which nothing is split, relative to size of the model.
    const double minLengthRatio = 1.0e-8;
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
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,"Node coordinates are not finite numbers.");
            }
        }
        for (const RTetraMeshData::Edge &rEdge : input.getEdges())
        {
            for (int nodeID : rEdge.nodeIDs)
            {
                if (nodeID < 0 || nodeID >= nNodes)
                {
                    throw RError(RError::Type::InvalidInput,R_ERROR_REF,"Edge refers to a node which does not exist (%d).",nodeID);
                }
            }
        }
        for (const RTetraMeshData::Triangle &rTriangle : input.getTriangles())
        {
            for (int nodeID : rTriangle.nodeIDs)
            {
                if (nodeID < 0 || nodeID >= nNodes)
                {
                    throw RError(RError::Type::InvalidInput,R_ERROR_REF,"Triangle refers to a node which does not exist (%d).",nodeID);
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
                        throw RError(RError::Type::InvalidInput,R_ERROR_REF,"Tetrahedron refers to a node which does not exist (%d).",nodeID);
                    }
                }
            }
        }
    }

    //! Find faces of input tetrahedra which lie on a boundary of a region.
    std::vector<std::array<int,3>> findRegionFaces(const RTetraMeshData &input, const std::vector<int> &nodeMap)
    {
        static const int faceNodes[4][3] = {{1,3,2},{0,2,3},{0,3,1},{0,1,2}};

        struct FaceInfo
        {
            std::array<int,3> nodeIDs;
            int marker;
            int count;
            bool boundary;
        };
        std::map<std::array<int,3>,FaceInfo> faces;

        for (const RTetraMeshData::Tetrahedron &rTetrahedron : input.getTetrahedra())
        {
            for (int i=0;i<4;i++)
            {
                const std::array<int,3> nodeIDs = {nodeMap[size_t(rTetrahedron.nodeIDs[faceNodes[i][0]])],
                                                   nodeMap[size_t(rTetrahedron.nodeIDs[faceNodes[i][1]])],
                                                   nodeMap[size_t(rTetrahedron.nodeIDs[faceNodes[i][2]])]};
                std::array<int,3> key = nodeIDs;
                std::sort(key.begin(),key.end());

                auto iter = faces.find(key);
                if (iter == faces.end())
                {
                    faces[key] = FaceInfo{nodeIDs,rTetrahedron.marker,1,false};
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

        std::vector<std::array<int,3>> regionFaces;
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
            for (int j=0;j<3 && !outside;j++)
            {
                double lower = std::numeric_limits<double>::max();
                double upper = std::numeric_limits<double>::lowest();
                for (int i=0;i<4;i++)
                {
                    const RTetraMeshData::Node &rNode = rNodes[size_t(rTetrahedron.nodeIDs[i])];
                    x[i][j] = j == 0 ? rNode.x : (j == 1 ? rNode.y : rNode.z);
                    lower = std::min(lower,x[i][j]);
                    upper = std::max(upper,x[i][j]);
                }
                outside = p[j] < lower || p[j] > upper;
            }
            if (outside)
            {
                continue;
            }
            const double orientation = RGeometricPredicates::orient3d(x[0],x[1],x[2],x[3]);
            if (orientation == 0.0)
            {
                continue;
            }
            const double sign = orientation > 0.0 ? 1.0 : -1.0;
            if (   sign * RGeometricPredicates::orient3d(x[1],x[3],x[2],p) >= 0.0
                && sign * RGeometricPredicates::orient3d(x[0],x[2],x[3],p) >= 0.0
                && sign * RGeometricPredicates::orient3d(x[0],x[3],x[1],p) >= 0.0
                && sign * RGeometricPredicates::orient3d(x[0],x[1],x[2],p) >= 0.0)
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
    generate(input,options,output,nullptr);
}

void RTetraMesher::generate(const RTetraMeshData &input, const RTetraMeshOptions &options, RTetraMeshData &output, RTetraMeshProfile *resultProfile)
{
    RTetraMeshProgress progress("Preparing volume mesh");
    RTetraMeshProfile profile;
    RTetraMeshProfile::Timer totalTimer(profile.total);
    const std::vector<RTetraMeshData::Node> &rInputNodes = input.getNodes();
    const int nInputNodes = int(rInputNodes.size());
    const bool reconstruct = options.getReconstruct() && !input.getTetrahedra().empty();

    validateInput(input,reconstruct);

    if (input.getTriangles().empty() && !reconstruct)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,"There are no triangles to generate tetrahedra from.");
    }

    // Delaunay tetrahedralization of input nodes.

    RTetraMesh mesh;
    mesh.setCoplanarTolerance(std::max(options.getTolerance(),1.0e-12));
    double lower[3] = {0.0,0.0,0.0};
    double upper[3] = {0.0,0.0,0.0};
    for (int i=0;i<nInputNodes;i++)
    {
        const double x[3] = {rInputNodes[size_t(i)].x,rInputNodes[size_t(i)].y,rInputNodes[size_t(i)].z};
        mesh.addPoint(x[0],x[1],x[2]);
        for (int j=0;j<3;j++)
        {
            lower[j] = i == 0 ? x[j] : std::min(lower[j],x[j]);
            upper[j] = i == 0 ? x[j] : std::max(upper[j],x[j]);
        }
    }
    const double modelSize = std::hypot(upper[0]-lower[0],upper[1]-lower[1],upper[2]-lower[2]);

    if (options.getVerbose())
    {
        RLogger::info("Tetrahedralizing %d nodes.\n",nInputNodes);
    }

    progress.initialize("Tetrahedralizing nodes");
    std::vector<int> duplicates;
    { RTetraMeshProfile::Timer timer(profile.tetrahedralization);
      RDelaunayTetrahedralizer::tetrahedralize(mesh,duplicates); }

    std::vector<int> nodeMap(static_cast<size_t>(nInputNodes));
    int nDuplicates = 0;
    for (int i=0;i<nInputNodes;i++)
    {
        nodeMap[size_t(i)] = duplicates[size_t(i)] >= 0 ? duplicates[size_t(i)] : i;
        if (duplicates[size_t(i)] >= 0)
        {
            nDuplicates++;
        }
    }
    if (nDuplicates > 0)
    {
        RLogger::warning("Found %d nodes with the same or nearly the same coordinates as other nodes, they are not used in generated mesh.\n",nDuplicates);
    }

    // Boundary.

    progress.initialize("Recovering mesh boundary");
    RTetraMeshProfile::Timer boundaryTimer(profile.boundary);
    RBoundaryRecovery boundary(mesh,options.getTolerance());
    boundary.setMinLength(minLengthRatio * modelSize);

    if (options.getUseSizeFunction())
    {
        for (int i=0;i<nInputNodes;i++)
        {
            const double size = rInputNodes[size_t(i)].size;
            if (size > 0.0 && boundary.getVertexSize(nodeMap[size_t(i)]) <= 0.0)
            {
                boundary.setVertexSize(nodeMap[size_t(i)],size);
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
        for (const std::array<int,3> &rFace : findRegionFaces(input,nodeMap))
        {
            boundary.addTriangle(rFace[0],rFace[1],rFace[2],RBoundaryRecovery::noMarker);
        }
    }
    for (const RTetraMeshData::Edge &rEdge : input.getEdges())
    {
        boundary.addEdge(nodeMap[size_t(rEdge.nodeIDs[0])],nodeMap[size_t(rEdge.nodeIDs[1])],rEdge.marker);
    }

    if (options.getVerbose())
    {
        RLogger::info("Recovering boundary.\n");
    }

    const uint64_t maxSteinerNodes = options.getMaxSteinerNodes() > 0 ? uint64_t(options.getMaxSteinerNodes()) : defaultMaxSteinerNodes;

    // Recovery which does not converge is stopped by a limit of its own, unless the limit given by options is lower.
    boundary.setMaxSteinerVertices(std::min(maxSteinerNodes,100 * uint64_t(nInputNodes) + 100000));
    boundary.initialize();
    boundary.recover();

    if (options.getVerbose())
    {
        RLogger::info("Boundary was recovered with %llu nodes added on edges and %llu nodes added on triangles.\n",
                      static_cast<unsigned long long>(boundary.getNSteinerVertices(RBoundaryRecovery::VertexSegment)),
                      static_cast<unsigned long long>(boundary.getNSteinerVertices(RBoundaryRecovery::VertexFacet)));
    }

    progress.initialize("Classifying mesh regions");
    int nRegions = boundary.classifyRegions();
    if (nRegions == 0)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,"Triangles do not enclose any volume.");
    }

    boundaryTimer.finish();
    RTetraMeshOptimizer optimizer(mesh,boundary,options);
    { RTetraMeshProfile::Timer timer(profile.initialRepair);optimizer.optimize(false); }

    // Refinement.

    if (options.getQualityMesh())
    {
        RTetraMeshProfile::Timer refinementTimer(profile.refinement);
        if (options.getVerbose())
        {
            RLogger::info("Refining mesh.\n");
        }

        // The limit is checked before a tetrahedron is refined, which may need more than one node.
        boundary.setMaxSteinerVertices(std::numeric_limits<uint64_t>::max());

        RTetraMeshRefiner refiner(mesh,boundary,options);
        refiner.setMinLength(minLengthRatio * modelSize);
        refiner.setMaxSteinerVertices(maxSteinerNodes);
        refiner.refine();
        profile.addWork(refiner.getProfile());
        if (options.getVerbose())
        {
            RLogger::info("Refinement safeguards rejected %llu interior insertions and %llu boundary splits; constrained cells may remain below requested quality.\n",
                          static_cast<unsigned long long>(refiner.getSpacingRejections()),static_cast<unsigned long long>(refiner.getBoundaryRejections()));
        }

        if (refiner.getStopped())
        {
            RLogger::warning("Mesh refinement was stopped after %llu nodes were added, some elements do not meet required quality.\n",
                             static_cast<unsigned long long>(boundary.getNSteinerVertices()));
        }

        progress.initialize("Recovering refined mesh boundary");
        boundary.recover();
        progress.initialize("Classifying refined mesh regions");
        nRegions = boundary.classifyRegions();
        if (nRegions == 0)
        {
            throw RError(RError::Type::Application,R_ERROR_REF,"Boundary was lost during mesh refinement.");
        }
    }

    if (options.getQualityMesh())
    {
        if (options.getVerbose()) RLogger::info("Optimizing mesh with bounded patch repairs.\n");
        RTetraMeshProfile::Timer timer(profile.optimization);
        optimizer.optimize();
    }

    // Output.

    progress.initialize("Preparing generated mesh output");
    output.clear();

    std::vector<RTetraMeshData::Node> &rNodes = output.getNodes();
    rNodes = rInputNodes;
    for (int i=nInputNodes;i<mesh.getNPoints();i++)
    {
        const double *p = mesh.getPoint(i);
        output.addNode(p[0],p[1],p[2],RTetraMeshData::noMarker,boundary.getVertexSize(i));
    }

    const int nSlots = mesh.getNTetraSlots();
    std::vector<int> tetraMap(static_cast<size_t>(nSlots),RTetraMeshData::noTetrahedron);
    int nTetrahedra = 0;
    for (int i=0;i<nSlots;i++)
    {
        if (mesh.isAlive(i) && !mesh.isGhost(i) && mesh.getTetra(i).region >= 0)
        {
            tetraMap[size_t(i)] = nTetrahedra++;
        }
    }

    // Markers of regions.
    std::vector<int> regionMarkers(static_cast<size_t>(nRegions));
    for (int i=0;i<nRegions;i++)
    {
        regionMarkers[size_t(i)] = i;
    }
    if (reconstruct)
    {
        int maxMarker = -1;
        for (const RTetraMeshData::Tetrahedron &rTetrahedron : input.getTetrahedra())
        {
            maxMarker = std::max(maxMarker,rTetrahedron.marker);
        }

        std::vector<std::map<int,int>> votes(static_cast<size_t>(nRegions));
        std::vector<int> nSamples(size_t(nRegions),0);
        for (int i=0;i<nSlots;i++)
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
            double center[3] = {0.0,0.0,0.0};
            for (int j=0;j<4;j++)
            {
                for (int k=0;k<3;k++)
                {
                    center[k] += 0.25 * mesh.getPoint(rTetra.v[j])[k];
                }
            }
            int marker;
            if (findInputMarker(input,center,marker))
            {
                votes[size_t(rTetra.region)][marker]++;
            }
        }
        for (int i=0;i<nRegions;i++)
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
    for (int i=0;i<nSlots;i++)
    {
        if (tetraMap[size_t(i)] == RTetraMeshData::noTetrahedron)
        {
            continue;
        }
        const RTetraMesh::Tetra &rTetra = mesh.getTetra(i);
        RTetraMeshData::Tetrahedron tetrahedron;
        for (int j=0;j<4;j++)
        {
            tetrahedron.nodeIDs[j] = rTetra.v[j];
            tetrahedron.neighborIDs[j] = tetraMap[size_t(rTetra.n[j])];
        }
        tetrahedron.marker = regionMarkers[size_t(rTetra.region)];
        rTetrahedra.push_back(tetrahedron);
    }

    for (const RBoundaryRecovery::SubFace &rFace : boundary.getSubFaces())
    {
        if (rFace.alive)
        {
            output.addTriangle(rFace.v[0],rFace.v[1],rFace.v[2],rFace.marker);
        }
    }
    for (const RBoundaryRecovery::SubSegment &rSegment : boundary.getSubSegments())
    {
        if (rSegment.alive && rSegment.marker != RBoundaryRecovery::noMarker)
        {
            output.addEdge(rSegment.v[0],rSegment.v[1],rSegment.marker);
        }
    }

    profile.addWork(optimizer.getProfile());
    progress.initialize("Checking mesh quality");
    RTetraMeshQuality::Summary quality;
    { RTetraMeshProfile::Timer timer(profile.diagnostics);quality=RTetraMeshQuality::summarize(output,options); }
    totalTimer.finish(); // Excludes optional report serialization and logging.
    if (!options.getQualityReportPath().empty())
    {
        progress.initialize("Writing mesh quality report");
        const RTetraMeshQuality::Work work={optimizer.getQueueVisits(),optimizer.getBudgetStops(),optimizer.getInsertedPoints(),optimizer.getOrderedVisits(),optimizer.getBoundaryPoints(),&profile};
        RTetraMeshQuality::writeReport(output,options,optimizer.getPatchAttempts(),options.getQualityReportPath(),&work);
    }
    if (resultProfile) *resultProfile=profile;
    if (options.getVerbose())
    {
        RLogger::info("Mesh phases (s): tetrahedralization %.3f, boundary %.3f, initial repair %.3f, refinement %.3f, optimization %.3f, diagnostics %.3f, total %.3f.\n",
                      profile.tetrahedralization,profile.boundary,profile.initialRepair,profile.refinement,profile.optimization,profile.diagnostics,profile.total);
        RLogger::info("Mesh work: %llu quality measurements, %llu quality cache hits, %llu star searches, %llu star cache hits, %llu local trial replacements.\n",
                      static_cast<unsigned long long>(profile.qualityMeasurements),static_cast<unsigned long long>(profile.qualityCacheHits),
                      static_cast<unsigned long long>(profile.starQueries),static_cast<unsigned long long>(profile.starCacheHits),static_cast<unsigned long long>(profile.trialReplacements));
        RLogger::info("Volume optimization visited %llu local cell patches, accepted %llu off-centers and reached %llu queue budgets.\n",
                      static_cast<unsigned long long>(optimizer.getQueueVisits()),static_cast<unsigned long long>(optimizer.getInsertedPoints()),static_cast<unsigned long long>(optimizer.getBudgetStops()));
        RLogger::info("Sub-degree angles: %llu elements have internal face pairs, %llu have prescribed face pairs (counts may overlap).\n",
                      static_cast<unsigned long long>(quality.freeBelowOne),static_cast<unsigned long long>(quality.prescribedBelowOne));
        RLogger::info("Mesh quality: minimum dihedral = %g degrees, radius/edge p95 = %g, maximum = %g, nearly flat = %llu.\n",
                      quality.minAngle,quality.radiusEdge95,quality.maxRadiusEdge,static_cast<unsigned long long>(quality.nearFlat));
    }
    if (options.getQualityMesh() && !quality.targetsMet())
    {
        RLogger::warning("Mesh quality targets remain unmet: %llu radius/edge, %llu dihedral, %llu volume/size, %llu nearly flat, %llu invalid elements (%llu elements have geometry-constrained angles).\n",
                         static_cast<unsigned long long>(quality.badRadius),static_cast<unsigned long long>(quality.badAngle),
                         static_cast<unsigned long long>(quality.tooLarge),static_cast<unsigned long long>(quality.nearFlat),
                         static_cast<unsigned long long>(quality.invalid),static_cast<unsigned long long>(quality.constrainedAngle));
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
    }
}
