#ifndef REL_TETRA_MESH_OPTIMIZER_H
#define REL_TETRA_MESH_OPTIMIZER_H
#include <array>
#include <vector>
#include <map>
#include <unordered_map>
#include "rel_tetra_mesh_quality.h"
#include "rel_tetra_mesh_profile.h"
#include "rel_boundary_recovery.h"
#include "rel_tetra_mesh_options.h"

//! Bounded, geometry-preserving volume flips and interior Steiner smoothing.
class RTetraMeshOptimizer
{
    RTetraMesh &mesh;
    RBoundaryRecovery &boundary;
    const RTetraMeshOptions &options;
    mutable RTetraMeshProfile profile;
    struct CachedElement { size_t stamp=0,boundaryEpoch=0;RTetraMeshQuality::Element value; };
    struct ArrayHash { size_t operator()(const std::array<int,4> &v) const noexcept { size_t h=0;for (int p : v) h^=size_t(p)+0x9e3779b9+(h<<6)+(h>>2);return h; } };
    mutable std::unordered_map<std::array<int,4>,CachedElement,ArrayHash> qualityCache;
    mutable std::vector<std::vector<int>> starCache;
    mutable std::vector<size_t> starCacheStamp;
    RTetraMeshQuality::Element element(const RTetraMesh &,const RBoundaryRecovery &,const int *vertices,int moved=-1,const double *position=nullptr) const;
    void findStar(int point,std::vector<int> &star) const;
    int findEdge(int a,int b) const;
    bool findFace(int a,int b,int c,RTetraMesh::Face &face) const;
    double objectiveScore(const RTetraMeshQuality::Element &quality) const;
    bool replace(const std::vector<int> &oldIDs, std::vector<std::array<int,4>> proposed, bool delaunay);
    bool flipFace(int a,int b,int c,bool delaunay);
    bool flipEdge(int a,int b,bool delaunay);
    bool relocate(int vertex);
    bool placePoint(int tetra);
    bool remeshPatch(int tetra);
    bool compoundFlip(int tetra);
    bool boundaryPatch(int tetra);
    std::vector<int> changed;
    std::map<std::array<int,4>,RTetraMeshQuality::PatchAttempts> attempts;
    RTetraMeshQuality::PatchAttempts rejected;
    size_t queueVisits = 0;
    size_t orderedVisits = 0;
    size_t budgetStops = 0;
    size_t insertedPoints = 0;
    size_t boundaryPoints = 0;
    size_t remeshAttempts = 0;
    size_t compoundAttempts = 0;
    size_t boundaryAttempts = 0;
    bool orderedRepairActive = false;
    bool simpleRepair = false;
    unsigned int orderedRepair(bool moveVertices = false);
    bool centroidRelocate(int vertex);
    size_t revision = 1;
    std::vector<size_t> vertexRevision;
    std::vector<size_t> geometryRevision;
    size_t boundaryEpoch=1;
    std::map<std::array<int,3>,size_t> edgeRevision;
    std::map<std::array<int,4>,size_t> faceRevision;
    void touch(const std::vector<int> &tetrahedra);
public:
    RTetraMeshOptimizer(RTetraMesh &mesh,RBoundaryRecovery &boundary,const RTetraMeshOptions &options);
    //! Both stages have finite pass and operation budgets. No boundary vertex is moved.
    unsigned int optimize(bool smoothVertices = true);
    //! Repair the neighborhood of freshly created tetrahedra with a bounded local sweep.
    unsigned int repairCreated(bool preferSimpleFlips = false);
    std::vector<RTetraMeshQuality::PatchAttempts> getPatchAttempts() const;
    const RTetraMeshProfile &getProfile() const { return profile; }
    size_t getQueueVisits() const { return queueVisits; }
    size_t getOrderedVisits() const { return orderedVisits; }
    size_t getBudgetStops() const { return budgetStops; }
    size_t getInsertedPoints() const { return insertedPoints; }
    size_t getBoundaryPoints() const { return boundaryPoints; }
};
#endif
