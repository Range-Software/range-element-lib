#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <queue>
#include <functional>

#include <rbl_error.h>

#include "rel_geometric_predicates.h"
#include "rel_delaunay_tetrahedralizer.h"
#include "rel_tetra_mesh_optimizer.h"
#include "rel_tetra_mesh_progress.h"
#include "rel_tetra_mesh_quality.h"

namespace
{

    double sizeViolation(const RBoundaryRecovery &boundary, const int *vertices, const RTetraMeshQuality::Element &q)
    {
        double size = 0;
        int count = 0;
        for (int i = 0; i < 4; i++)
        {
            if (boundary.getVertexSize(vertices[i]) > 0)
            {
                size += boundary.getVertexSize(vertices[i]);
                count++;
            }
        }
        return count > 0 ? q.circumradius / (.75 * size / count) : 0.0;
    }

    RTetraMeshQuality::Element uncachedElement(const RTetraMesh &mesh, const RBoundaryRecovery &boundary, const int *v, int moved = -1, const double *position = nullptr)
    {
        const double *p[4];
        bool fixed[4];
        for (int i = 0; i < 4; i++)
        {
            p[i] = v[i] == moved ? position : mesh.getPoint(v[i]);
            const int *f = RTetraMesh::faceVertices[i];
            fixed[i] = boundary.findSubFace(v[f[0]], v[f[1]], v[f[2]]) >= 0;
        }
        return RTetraMeshQuality::measure(p[0], p[1], p[2], p[3], fixed);
    }

    //! Sort and remove repeated keys: the order and content of an ordered set, without its nodes.
    template <typename T>
    void makeUnique(std::vector<T> &keys)
    {
        std::sort(keys.begin(), keys.end());
        keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
    }

    inline bool hasVertex(const int *v, int vertex)
    {
        return v[0] == vertex || v[1] == vertex || v[2] == vertex || v[3] == vertex;
    }

    double qualityScore(const RTetraMeshQuality::Element &q)
    {
        return std::min({q.meanRatio, 1.0 / (1.0 + q.radiusEdge * q.radiusEdge), q.minFreeAngle / 180.0});
    }

}

RTetraMeshOptimizer::RTetraMeshOptimizer(RTetraMesh &mesh, RBoundaryRecovery &boundary, const RTetraMeshOptions &options)
    : mesh(mesh)
    , boundary(boundary)
    , options(options)
    , vertexRevision(size_t(mesh.getNPoints()), 1)
    , geometryRevision(size_t(mesh.getNPoints()), 1)
{
}

RTetraMeshQuality::Element RTetraMeshOptimizer::element(const RTetraMesh &m, const RBoundaryRecovery &b, const int *v, int moved, const double *position) const
{
    // Cache only the live mesh, retaining vertex order to preserve rounding and
    // predicate behavior. Trial positions and independent boundary copies bypass it.
    if (&m != &mesh || &b != &boundary || moved >= 0)
    {
        profile.qualityMeasurements++;
        return uncachedElement(m, b, v, moved, position);
    }
    size_t stamp = 1;
    for (int i = 0; i < 4; i++)
    {
        if (size_t(v[i]) >= geometryRevision.size())
        {
            profile.qualityMeasurements++;
            return uncachedElement(m, b, v);
        }
        stamp = std::max(stamp, geometryRevision[size_t(v[i])]);
    }
    const size_t cacheSize = 65536;
    if (qualityCache.empty())
    {
        qualityCache.resize(cacheSize);
    }
    const std::array<int, 4> key = {v[0], v[1], v[2], v[3]};
    auto &entry = qualityCache[KeyHash()(key) & (cacheSize - 1)];
    if (entry.key == key && entry.stamp == stamp && entry.boundaryEpoch == boundaryEpoch)
    {
        profile.qualityCacheHits++;
        return entry.value;
    }
    profile.qualityMeasurements++;
    entry.key = key;
    entry.stamp = stamp;
    entry.boundaryEpoch = boundaryEpoch;
    return entry.value = uncachedElement(m, b, v);
}

void RTetraMeshOptimizer::findStar(int point, std::vector<int> &star) const
{
    starCache.resize(size_t(mesh.getNPoints()));
    starCacheStamp.resize(size_t(mesh.getNPoints()), 0);
    if (starCacheStamp[size_t(point)] == vertexRevision[size_t(point)])
    {
        profile.starCacheHits++;
        star = starCache[size_t(point)];
        return;
    }
    profile.starQueries++;
    mesh.findStar(point, star);
    starCache[size_t(point)] = star;
    starCacheStamp[size_t(point)] = vertexRevision[size_t(point)];
}

const std::vector<int> &RTetraMeshOptimizer::starOf(int point) const
{
    starCache.resize(size_t(mesh.getNPoints()));
    starCacheStamp.resize(size_t(mesh.getNPoints()), 0);
    auto &star = starCache[size_t(point)];
    if (starCacheStamp[size_t(point)] == vertexRevision[size_t(point)])
    {
        profile.starCacheHits++;
        return star;
    }
    profile.starQueries++;
    mesh.findStar(point, star);
    starCacheStamp[size_t(point)] = vertexRevision[size_t(point)];
    return star;
}

int RTetraMeshOptimizer::findEdge(int a, int b) const
{
    for (int id : starOf(a))
    {
        if (hasVertex(mesh.getTetra(id).v, b))
        {
            return id;
        }
    }
    return RTetraMesh::noTetra;
}

bool RTetraMeshOptimizer::findFace(int a, int b, int c, RTetraMesh::Face &face) const
{
    for (int id : starOf(a))
    {
        const auto &t = mesh.getTetra(id);
        int count = 0, opposite = -1;
        for (int i = 0; i < 4; i++)
        {
            if (t.v[i] == a || t.v[i] == b || t.v[i] == c)
            {
                count++;
            }
            else
            {
                opposite = i;
            }
        }
        if (count == 3)
        {
            face = {id, opposite};
            return true;
        }
    }
    return false;
}

double RTetraMeshOptimizer::objectiveScore(const RTetraMeshQuality::Element &quality) const
{
    auto q = quality;
    if (orderedRepairActive || simpleRepair)
    {
        q.minFreeAngle = q.minAngle;
    }
    return qualityScore(q);
}

struct RTetraMeshOptimizer::PatchState
{
    enum Failure
    {
        None,
        Region,
        Degenerate,
        Topology
    };

    bool ready = false;
    Failure failure = None;
    double center[3] = {0, 0, 0}, scale = 0;
    double score = 0, energy = 0, volume = 0, size = 0, angle = 0, ratio = 0;
    size_t subOne = 0;
};

bool RTetraMeshOptimizer::replace(const std::vector<int> &oldIDs, std::vector<std::array<int, 4>> &proposed, bool delaunay, bool boundaryChecked, PatchState *shared)
{
    if (oldIDs.empty())
    {
        return false;
    }
    const int region = this->mesh.getTetra(oldIDs.front()).region;
    std::set<std::array<int, 3>> newFaces;
    std::set<std::array<int, 2>> newEdges;
    for (auto &v : proposed)
    {
        double sign = RGeometricPredicates::orient3d(mesh.getPoint(v[0]), mesh.getPoint(v[1]), mesh.getPoint(v[2]), mesh.getPoint(v[3]));
        if (sign < 0.0)
        {
            std::swap(v[0], v[1]);
        }
        if (sign == 0.0)
        {
            rejected.topology++;
            return false;
        }
        if (boundaryChecked)
        {
            continue;
        }
        for (int i = 0; i < 4; i++)
        {
            const int *f = RTetraMesh::faceVertices[i];
            std::array<int, 3> key = {v[f[0]], v[f[1]], v[f[2]]};
            std::sort(key.begin(), key.end());
            newFaces.insert(key);
            for (int j = i + 1; j < 4; j++)
            {
                newEdges.insert({std::min(v[i], v[j]), std::max(v[i], v[j])});
            }
        }
    }
    // Everything measured on the old cells is independent of the proposal once the
    // boundary is known to be preserved, so one patch state serves all its proposals.
    PatchState local;
    PatchState &state = shared && boundaryChecked ? *shared : local;
    auto fail = [&](PatchState::Failure failure)
    {
        state.failure = failure;
        if (failure == PatchState::Region)
        {
            rejected.region++;
        }
        if (failure == PatchState::Topology)
        {
            rejected.topology++;
        }
        return false;
    };
    auto metrics = [&](const std::vector<std::array<int, 4>> &tets, double &worst, double &energy, double &maxVolume, double &sizeExcess, double &angle, double &ratio, size_t &subOne, double scoreFloor = -1.0)
    {
        angle = 180;
        ratio = 0;
        subOne = 0;
        worst = 1.0;
        energy = 0;
        maxVolume = 0;
        sizeExcess = 0;
        for (const auto &v : tets)
        {
            const auto q = element(mesh, boundary, v.data());
            if (!(q.volume > 0.0) || !std::isfinite(q.volume))
            {
                return false;
            }
            angle = std::min(angle, q.minFreeAngle);
            ratio = std::max(ratio, q.radiusEdge);
            subOne += q.minFreeAngle < 1.0;
            worst = std::min(worst, objectiveScore(q));
            // A bottleneck below this bound cannot recover in later cells.
            if (worst <= scoreFloor)
            {
                return false;
            }
            maxVolume = std::max(maxVolume, q.volume);
            if (options.getUseSizeFunction())
            {
                sizeExcess = std::max(sizeExcess, sizeViolation(boundary, v.data(), q));
            }
            if (delaunay)
            {
                double sum = 0;
                for (int id : v)
                {
                    for (int i = 0; i < 3; i++)
                    {
                        sum += std::pow((mesh.getPoint(id)[i] - state.center[i]) / state.scale, 2);
                    }
                }
                energy += (q.volume / (state.scale * state.scale * state.scale)) * sum;
            }
        }
        return true;
    };
    if (state.ready)
    {
        if (state.failure != PatchState::None)
        {
            return fail(state.failure);
        }
    }
    else
    {
        state.ready = true;
        std::vector<int> points;
        points.reserve(4 * oldIDs.size());
        for (int id : oldIDs)
        {
            if (!mesh.isAlive(id) || mesh.isGhost(id) || mesh.getTetra(id).region != region)
            {
                return fail(PatchState::Region);
            }
            const auto &t = mesh.getTetra(id);
            points.insert(points.end(), t.v, t.v + 4);
            if (boundaryChecked)
            {
                continue;
            }
            for (int i = 0; i < 4; i++)
            {
                const int *f = RTetraMesh::faceVertices[i];
                std::array<int, 3> key = {t.v[f[0]], t.v[f[1]], t.v[f[2]]};
                std::sort(key.begin(), key.end());
                if (boundary.findSubFace(key[0], key[1], key[2]) >= 0 && !newFaces.count(key))
                {
                    rejected.boundary++;
                    return false;
                }
                for (int j = i + 1; j < 4; j++)
                {
                    if (boundary.findSubSegment(t.v[i], t.v[j]) >= 0 && !newEdges.count({std::min(t.v[i], t.v[j]), std::max(t.v[i], t.v[j])}))
                    {
                        rejected.boundary++;
                        return false;
                    }
                }
            }
        }
        // Ascending order keeps the rounding of the center independent of cell order.
        std::sort(points.begin(), points.end());
        points.erase(std::unique(points.begin(), points.end()), points.end());
        for (int v : points)
        {
            for (int i = 0; i < 3; i++)
            {
                state.center[i] += mesh.getPoint(v)[i] / double(points.size());
            }
        }
        for (int v : points)
        {
            for (int i = 0; i < 3; i++)
            {
                state.scale = std::max(state.scale, std::fabs(mesh.getPoint(v)[i] - state.center[i]));
            }
        }
        if (!(state.scale > 0.0))
        {
            return fail(PatchState::Degenerate);
        }
        std::vector<std::array<int, 4>> old;
        old.reserve(oldIDs.size());
        for (int id : oldIDs)
        {
            const auto &t = mesh.getTetra(id);
            old.push_back({t.v[0], t.v[1], t.v[2], t.v[3]});
        }
        if (!metrics(old, state.score, state.energy, state.volume, state.size, state.angle, state.ratio, state.subOne))
        {
            return fail(PatchState::Topology);
        }
    }
    const double before = state.score, energyBefore = state.energy, volumeBefore = state.volume, sizeBefore = state.size, angleBefore = state.angle, ratioBefore = state.ratio;
    const size_t subBefore = state.subOne;
    double after, energyAfter, volumeAfter, sizeAfter, angleAfter, ratioAfter;
    size_t subAfter;
    const double scoreFloor = orderedRepairActive && !delaunay ? before * (1 + 1e-5) + 1e-12 : -1.0;
    if (!metrics(proposed, after, energyAfter, volumeAfter, sizeAfter, angleAfter, ratioAfter, subAfter, scoreFloor))
    {
        if (scoreFloor >= 0)
        {
            rejected.quality++;
        }
        else
        {
            rejected.topology++;
        }
        return false;
    }
    if (!delaunay && !orderedRepairActive && subAfter > subBefore)
    {
        rejected.quality++;
        return false;
    }
    if (!delaunay && !orderedRepairActive && (angleAfter + 1e-8 < std::min(angleBefore, 1.0) || ratioAfter > std::max(ratioBefore, options.getRadiusEdgeRatio()) * (1.0 + 1e-8)))
    {
        rejected.quality++;
        return false;
    }
    if (sizeAfter > std::max(1.0, sizeBefore) * (1.0 + 1e-10))
    {
        return false;
    }
    if (options.getVolumeConstraint() > 0.0 && volumeAfter > std::max(volumeBefore, options.getVolumeConstraint()) * (1.0 + 1e-10))
    {
        return false;
    }
    if (delaunay)
    {
        if (!(energyAfter < energyBefore * (1.0 - 1e-11)) || after < std::max(1e-10, before * .25))
        {
            rejected.quality++;
            return false;
        }
    }
    else
    {
        if (orderedRepairActive)
        {
            if (!(after > before * (1 + 1e-5) + 1e-12))
            {
                rejected.quality++;
                return false;
            }
        }
        else
        {
            const double target = std::max(1.0, options.getMinDihedralAngle());
            const double oldAngle = std::min(angleBefore, target), newAngle = std::min(angleAfter, target);
            const bool betterAngle = newAngle > oldAngle + 1e-7;
            if (newAngle + 1e-7 < oldAngle || (!betterAngle && !(after > before * (1.0 + 1e-5) + 1e-12)) || (betterAngle && after < before * .75))
            {
                rejected.quality++;
                return false;
            }
        }
        // Improve the clipped repairable-angle objective first, then shape.
        // Angular improvements may spend at most 25% of the local shape score;
        // explicit radius/angle safety floors and finite budgets remain binding.
    }
    if (!mesh.replaceTetras(oldIDs, proposed))
    {
        rejected.topology++;
        return false;
    }
    changed = mesh.getCreatedTetras();
    profile.acceptedChanges++;
    touch(changed);
    auto &created = boundary.getCreatedTetras();
    created.insert(created.end(), mesh.getCreatedTetras().begin(), mesh.getCreatedTetras().end());
    return true;
}

bool RTetraMeshOptimizer::flipFace(int a, int b, int c, bool delaunay)
{
    RTetraMesh::Face face;
    if (boundary.findSubFace(a, b, c) >= 0)
    {
        rejected.boundary++;
        return false;
    }
    if (!findFace(a, b, c, face))
    {
        return false;
    }
    const auto t = mesh.getTetra(face.tetra);
    const int other = t.n[face.face];
    if (mesh.isGhost(face.tetra) || mesh.isGhost(other))
    {
        return false;
    }
    size_t stamp = 1;
    for (int v : t.v)
    {
        stamp = std::max(stamp, vertexRevision[size_t(v)]);
    }
    for (int v : mesh.getTetra(other).v)
    {
        stamp = std::max(stamp, vertexRevision[size_t(v)]);
    }
    const std::array<int, 4> cacheKey = {a, b, c, int(delaunay) | (int(orderedRepairActive) << 1) | (int(simpleRepair) << 2)};
    if (faceRevision[cacheKey] == stamp)
    {
        profile.queueSkips++;
        return false;
    }
    faceRevision[cacheKey] = stamp;
    int e = -1;
    for (int v : mesh.getTetra(other).v)
    {
        if (v != a && v != b && v != c)
        {
            e = v;
        }
    }
    const int d = t.v[face.face];
    if (e < 0)
    {
        return false;
    }
    // The Delaunay objective is the volume under the lifted mesh, which a 2-3 flip
    // lowers only if the opposite vertex lies inside the circumscribed sphere.
    if (delaunay && !(RGeometricPredicates::inSphere(mesh.getPoint(t.v[0]), mesh.getPoint(t.v[1]), mesh.getPoint(t.v[2]), mesh.getPoint(t.v[3]), mesh.getPoint(e)) > 0.0))
    {
        profile.delaunaySkips++;
        return false;
    }
    if (findEdge(d, e) != RTetraMesh::noTetra)
    {
        return false;
    }
    // The only prescribed entity a 2-3 flip can remove is the face tested above.
    std::vector<std::array<int, 4>> proposed = {{d, e, a, b}, {d, e, b, c}, {d, e, c, a}};
    return replace({face.tetra, other}, proposed, delaunay, true);
}

bool RTetraMeshOptimizer::flipEdge(int a, int b, bool delaunay)
{
    if (boundary.findSubSegment(a, b) >= 0)
    {
        rejected.boundary++;
        return false;
    }
    std::vector<int> ids;
    for (int id : starOf(a))
    {
        const auto &t = mesh.getTetra(id);
        if (!hasVertex(t.v, b))
        {
            continue;
        }
        if (mesh.isGhost(id))
        {
            return false;
        }
        ids.push_back(id);
    }
    if (ids.size() < 3 || ids.size() > ((orderedRepairActive || simpleRepair) ? 8 : 10))
    {
        return false;
    }

    // Vertices around the edge with their neighbors in the order cells are met.
    // At most ten cells remain here, so a linear search beats a tree.
    struct RingVertex
    {
        int vertex, next[2], count;
    };

    std::array<RingVertex, 20> ring;
    size_t ringSize = 0;
    auto link = [&](int x, int y)
    {
        for (size_t i = 0; i < ringSize; i++)
        {
            if (ring[i].vertex == x)
            {
                if (ring[i].count < 2)
                {
                    ring[i].next[ring[i].count] = y;
                }
                ring[i].count++;
                return;
            }
        }
        ring[ringSize++] = {x, {y, -1}, 1};
    };
    for (int id : ids)
    {
        const auto &t = mesh.getTetra(id);
        int x = -1, y = -1;
        for (int v : t.v)
        {
            if (v != a && v != b)
            {
                if (x < 0)
                {
                    x = v;
                }
                else
                {
                    y = v;
                }
            }
        }
        link(x, y);
        link(y, x);
    }
    size_t stamp = 1;
    vertexRevision.resize(size_t(mesh.getNPoints()), revision);
    for (int id : ids)
    {
        for (int v : mesh.getTetra(id).v)
        {
            stamp = std::max(stamp, vertexRevision[size_t(v)]);
        }
    }
    const std::array<int, 3> cacheKey = {a, b, int(delaunay) | (int(orderedRepairActive) << 1) | (int(simpleRepair) << 2)};
    if (edgeRevision[cacheKey] == stamp)
    {
        profile.queueSkips++;
        return false;
    }
    edgeRevision[cacheKey] = stamp;
    if (delaunay)
    {
        // The lifted mesh is convex across every locally Delaunay face. If all faces
        // at the edge are, it lies below any other mesh of these cells' vertices
        // and no removal of the edge can lower the Delaunay objective.
        bool improvable = false;
        for (size_t i = 0; i < ids.size() && !improvable; i++)
        {
            const auto &t = mesh.getTetra(ids[i]);
            for (int j = 0; j < 4 && !improvable; j++)
            {
                // The face opposite to a vertex of the ring is shared with the next cell at the edge.
                if (t.v[j] == a || t.v[j] == b || t.n[j] < ids[i])
                {
                    continue;
                }
                const auto &n = mesh.getTetra(t.n[j]);
                for (int v : n.v)
                {
                    if (v == RTetraMesh::ghostVertex || hasVertex(t.v, v))
                    {
                        continue;
                    }
                    improvable = RGeometricPredicates::inSphere(mesh.getPoint(t.v[0]), mesh.getPoint(t.v[1]), mesh.getPoint(t.v[2]), mesh.getPoint(t.v[3]), mesh.getPoint(v)) > 0.0;
                }
            }
        }
        if (!improvable)
        {
            profile.delaunaySkips++;
            return false;
        }
    }
    if (ringSize != ids.size())
    {
        return false;
    }
    for (size_t i = 0; i < ringSize; i++)
    {
        if (ring[i].count != 2)
        {
            return false;
        }
    }
    std::vector<int> cycle;
    cycle.reserve(ringSize);
    // The cycle starts at the lowest vertex, as it did when the ring was an ordered map.
    int previous = -1, current = ring[0].vertex;
    for (size_t i = 1; i < ringSize; i++)
    {
        current = std::min(current, ring[i].vertex);
    }
    for (size_t i = 0; i < ringSize; i++)
    {
        if (std::find(cycle.begin(), cycle.end(), current) != cycle.end())
        {
            return false;
        }
        cycle.push_back(current);
        const int *next = nullptr;
        for (size_t j = 0; j < ringSize; j++)
        {
            if (ring[j].vertex == current)
            {
                next = ring[j].next;
            }
        }
        int n = next[0] == previous ? next[1] : next[0];
        previous = current;
        current = n;
    }
    if (current != cycle.front())
    {
        return false;
    }
    for (int v : cycle)
    {
        if (boundary.findSubFace(a, b, v) >= 0)
        {
            rejected.boundary++;
            return false;
        }
    }
    // Removing the edge keeps every other face and edge of its cells, so the
    // boundary needs no further test and the old cells are measured only once.
    PatchState patch;
    auto fans = [&]()
    {
        for (size_t root = 0; root < cycle.size(); root++)
        {
            std::vector<std::array<int, 4>> proposed;
            bool allowed = true;
            for (size_t j = 1; j + 1 < cycle.size(); j++)
            {
                const int x = cycle[root], y = cycle[(root + j) % cycle.size()], z = cycle[(root + j + 1) % cycle.size()];
                if (j > 1 && findEdge(x, y) != RTetraMesh::noTetra)
                {
                    allowed = false;
                    break;
                }
                proposed.push_back({x, y, z, a});
                proposed.push_back({x, z, y, b});
            }
            if (allowed && replace(ids, proposed, delaunay, true, &patch))
            {
                return true;
            }
        }
        return false;
    };
    if (orderedRepairActive || simpleRepair)
    {
        return fans();
    }
    // Bottleneck dynamic programming considers all polygon triangulations,
    // including non-fan triangulations, in O(n^3) time and O(n^2) storage.
    const int n = int(cycle.size());
    double oldScore = 1, oldAngle = 180, oldRatio = 0, oldVolume = 0, oldSize = 0;
    double centroid[3] = {0, 0, 0};
    std::set<int> patchPoints = {a, b};
    patchPoints.insert(cycle.begin(), cycle.end());
    for (int v : patchPoints)
    {
        for (int i = 0; i < 3; i++)
        {
            centroid[i] += mesh.getPoint(v)[i] / double(patchPoints.size());
        }
    }
    for (int id : ids)
    {
        const auto &t = mesh.getTetra(id);
        const auto q = element(mesh, boundary, t.v);
        oldScore = std::min(oldScore, objectiveScore(q));
        oldAngle = std::min(oldAngle, q.minFreeAngle);
        oldRatio = std::max(oldRatio, q.radiusEdge);
        oldVolume = std::max(oldVolume, q.volume);
        if (options.getUseSizeFunction())
        {
            oldSize = std::max(oldSize, sizeViolation(boundary, t.v, q));
        }
    }
    auto parabolic = [&](const int *v, const RTetraMeshQuality::Element &q)
    {
        double sum = 0;
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                sum += std::pow(mesh.getPoint(v[i])[j] - centroid[j], 2);
            }
        }
        return q.volume * sum;
    };
    double best[10][10] = {}, total[10][10] = {}, energy[10][10] = {}, angles[10][10] = {};
    int split[10][10];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            best[i][j] = -1;
            split[i][j] = -1;
        }
    }
    for (int i = 0; i + 1 < n; i++)
    {
        best[i][i + 1] = 1;
        angles[i][i + 1] = 180;
    }
    auto diagonalAllowed = [&](int i, int j)
    {
        return j == i + 1 || (i == 0 && j == n - 1) || findEdge(cycle[i], cycle[j]) == RTetraMesh::noTetra;
    };
    for (int width = 2; width < n; width++)
    {
        for (int i = 0; i + width < n; i++)
        {
            const int j = i + width;
            if (!diagonalAllowed(i, j))
            {
                continue;
            }
            for (int k = i + 1; k < j; k++)
            {
                if (best[i][k] < 0 || best[k][j] < 0)
                {
                    continue;
                }
                const int x = cycle[i], y = cycle[k], z = cycle[j];
                double sa = RGeometricPredicates::orient3d(mesh.getPoint(x), mesh.getPoint(y), mesh.getPoint(z), mesh.getPoint(a));
                double sb = RGeometricPredicates::orient3d(mesh.getPoint(x), mesh.getPoint(y), mesh.getPoint(z), mesh.getPoint(b));
                if (sa == 0 || sb == 0 || (sa > 0) == (sb > 0))
                {
                    continue;
                }
                const int va[4] = {x, sa > 0 ? y : z, sa > 0 ? z : y, a}, vb[4] = {x, sb > 0 ? y : z, sb > 0 ? z : y, b};
                const auto qa = element(mesh, boundary, va), qb = element(mesh, boundary, vb);
                if (!delaunay && !orderedRepairActive && (std::min(qa.minFreeAngle, qb.minFreeAngle) + 1e-8 < std::min(oldAngle, 1.0) || std::max(qa.radiusEdge, qb.radiusEdge) > std::max(oldRatio, options.getRadiusEdgeRatio()) * (1 + 1e-8)))
                {
                    continue;
                }
                if (options.getVolumeConstraint() > 0 && std::max(qa.volume, qb.volume) > std::max(oldVolume, options.getVolumeConstraint()) * (1 + 1e-10))
                {
                    continue;
                }
                if (options.getUseSizeFunction() && std::max(sizeViolation(boundary, va, qa), sizeViolation(boundary, vb, qb)) > std::max(1.0, oldSize) * (1 + 1e-10))
                {
                    continue;
                }
                if (delaunay && std::min(objectiveScore(qa), objectiveScore(qb)) < std::max(1e-10, oldScore * .25))
                {
                    continue;
                }
                const double score = std::min({best[i][k], best[k][j], objectiveScore(qa), objectiveScore(qb)});
                const double sum = total[i][k] + total[k][j] + objectiveScore(qa) + objectiveScore(qb);
                const double cost = energy[i][k] + energy[k][j] + parabolic(va, qa) + parabolic(vb, qb);
                const double angle = std::min({angles[i][k], angles[k][j], qa.minFreeAngle, qb.minFreeAngle, std::max(1.0, options.getMinDihedralAngle())});
                const bool scoreImproves = score > best[i][j] + 1e-12 || (std::fabs(score - best[i][j]) <= 1e-12 && sum > total[i][j]);
                const bool shapeImproves = delaunay ? scoreImproves : (angle > angles[i][j] + 1e-7 || (std::fabs(angle - angles[i][j]) <= 1e-7 && scoreImproves));
                const bool sameEnergy = std::fabs(cost - energy[i][j]) <= 1e-12 * std::max(std::fabs(cost), std::fabs(energy[i][j]));
                const bool improves = delaunay ? (cost < energy[i][j] * (1 - 1e-12) || (sameEnergy && shapeImproves)) : shapeImproves;
                if (split[i][j] < 0 || improves)
                {
                    best[i][j] = score;
                    total[i][j] = sum;
                    energy[i][j] = cost;
                    angles[i][j] = angle;
                    split[i][j] = k;
                }
            }
        }
    }
    if (split[0][n - 1] < 0)
    {
        rejected.topology++;
    }
    std::vector<std::array<int, 4>> proposed;
    std::function<void(int, int)> triangulate = [&](int i, int j)
    {
        if (j <= i + 1)
        {
            return;
        }
        const int k = split[i][j];
        proposed.push_back({cycle[i], cycle[k], cycle[j], a});
        proposed.push_back({cycle[i], cycle[j], cycle[k], b});
        triangulate(i, k);
        triangulate(k, j);
    };
    if (split[0][n - 1] >= 0)
    {
        triangulate(0, n - 1);
        if (replace(ids, proposed, delaunay, true, &patch))
        {
            return true;
        }
    }
    // The bottleneck solution can be geometrically inadmissible for a nonplanar
    // link, or conflict with an existing external face. Retain bounded alternative
    // fans rather than letting that one rejected proposal hide all legal repairs.
    for (int root = 0; root < n; root++)
    {
        proposed.clear();
        bool allowed = true;
        for (int j = 1; j + 1 < n; j++)
        {
            const int x = cycle[size_t(root)], y = cycle[size_t((root + j) % n)], z = cycle[size_t((root + j + 1) % n)];
            if (j > 1 && findEdge(x, y) != RTetraMesh::noTetra)
            {
                allowed = false;
                break;
            }
            proposed.push_back({x, y, z, a});
            proposed.push_back({x, z, y, b});
        }
        if (allowed && replace(ids, proposed, delaunay, true, &patch))
        {
            return true;
        }
    }
    return false;
}

void RTetraMeshOptimizer::touch(const std::vector<int> &ids)
{
    vertexRevision.resize(size_t(mesh.getNPoints()), revision);
    geometryRevision.resize(size_t(mesh.getNPoints()), revision);
    ++revision;
    for (int id : ids)
    {
        if (mesh.isAlive(id))
        {
            for (int v : mesh.getTetra(id).v)
            {
                if (v >= 0)
                {
                    vertexRevision[size_t(v)] = revision;
                }
            }
        }
    }
}

bool RTetraMeshOptimizer::relocate(int id)
{
    RTetraMeshProfile::Timer timer(profile.relocation);
    if (boundary.getVertexType(id) != RBoundaryRecovery::VertexVolume || boundary.isBoundaryVertex(id))
    {
        return false;
    }
    std::vector<int> star;
    findStar(id, star);
    if (star.empty())
    {
        return false;
    }
    std::set<int> neighbors;
    double old[3], center[3] = {0, 0, 0};
    std::copy(mesh.getPoint(id), mesh.getPoint(id) + 3, old);
    double before = 1, maxVolume = 0, sizeBefore = 0, angleBefore = 180, ratioBefore = 0, scale = std::numeric_limits<double>::max();
    for (int t : star)
    {
        if (mesh.isGhost(t) || mesh.getTetra(t).region < 0)
        {
            return false;
        }
        const auto &v = mesh.getTetra(t);
        for (int point : v.v)
        {
            if (point != id)
            {
                neighbors.insert(point);
            }
        }
        const auto q = element(mesh, boundary, v.v);
        before = std::min(before, objectiveScore(q));
        maxVolume = std::max(maxVolume, q.volume);
        angleBefore = std::min(angleBefore, q.minFreeAngle);
        ratioBefore = std::max(ratioBefore, q.radiusEdge);
        if (options.getUseSizeFunction())
        {
            sizeBefore = std::max(sizeBefore, sizeViolation(boundary, v.v, q));
        }
    }
    for (int point : neighbors)
    {
        double length = 0;
        for (int i = 0; i < 3; i++)
        {
            center[i] += mesh.getPoint(point)[i] / double(neighbors.size());
            length += std::pow(mesh.getPoint(point)[i] - old[i], 2);
        }
        scale = std::min(scale, std::sqrt(length));
    }

    struct StarCell
    {
        std::array<int, 4> vertices;
        std::array<const double *, 4> points;
        std::array<bool, 4> fixed;
        RTetraMeshQuality::Element previous;
        int moving;
    };

    std::vector<StarCell> cells;
    cells.reserve(star.size());
    for (int t : star)
    {
        const auto &v = mesh.getTetra(t);
        StarCell cell;
        cell.previous = element(mesh, boundary, v.v);
        for (int i = 0; i < 4; i++)
        {
            cell.vertices[i] = v.v[i];
            cell.points[i] = mesh.getPoint(v.v[i]);
            if (v.v[i] == id)
            {
                cell.moving = i;
            }
            const int *f = RTetraMesh::faceVertices[i];
            cell.fixed[i] = boundary.findSubFace(v.v[f[0]], v.v[f[1]], v.v[f[2]]) >= 0;
        }
        cells.push_back(cell);
    }
    auto objective = [&](const double *p)
    {
        double worst = 1, angle = 180;
        for (const auto &cell : cells)
        {
            auto points = cell.points;
            points[size_t(cell.moving)] = p;
            // Reject invalid orientation before expensive radius and angle evaluation.
            if (!(RGeometricPredicates::orient3d(points[0], points[1], points[2], points[3]) > 0))
            {
                return std::array<double, 2>{-1, -1};
            }
            profile.qualityMeasurements++;
            const auto quality = RTetraMeshQuality::measure(points[0], points[1], points[2], points[3], cell.fixed.data());
            const auto &previous = cell.previous;
            if (quality.minFreeAngle + 1e-8 < std::min(previous.minFreeAngle, 1.0))
            {
                return std::array<double, 2>{-1, -1};
            }
            if (!(quality.volume > 0) || quality.minFreeAngle + 1e-8 < std::min(angleBefore, 1.0) || quality.radiusEdge > std::max(ratioBefore, options.getRadiusEdgeRatio()) * (1 + 1e-8))
            {
                return std::array<double, 2>{-1, -1};
            }
            if (options.getVolumeConstraint() > 0 && quality.volume > std::max(maxVolume, options.getVolumeConstraint()) * (1 + 1e-10))
            {
                return std::array<double, 2>{-1, -1};
            }
            if (options.getUseSizeFunction() && sizeViolation(boundary, cell.vertices.data(), quality) > std::max(1.0, sizeBefore) * (1 + 1e-10))
            {
                return std::array<double, 2>{-1, -1};
            }
            worst = std::min(worst, objectiveScore(quality));
            angle = std::min(angle, quality.minFreeAngle);
        }
        return std::array<double, 2>{std::min(angle, std::max(1.0, options.getMinDihedralAngle())), worst};
    };
    double best[3];
    std::copy(old, old + 3, best);
    auto bestValue = std::array<double, 2>{std::min(angleBefore, std::max(1.0, options.getMinDihedralAngle())), before};
    auto improves = [](const auto &a, const auto &b)
    {
        return (a[0] > b[0] + 1e-7 && a[1] >= b[1] * .75) || (std::fabs(a[0] - b[0]) <= 1e-7 && a[1] > b[1] * (1 + 1e-5) + 1e-12);
    };
    // Bounded nonsmooth trust-region search optimizes the actual worst element,
    // rather than assuming that the neighbor centroid is an improving direction.
    for (double radius : {.25, .1, .04})
    {
        double origin[3];
        std::copy(best, best + 3, origin);
        for (int direction = 0; direction < (angleBefore < 1.0 ? 29 : 9); direction++)
        {
            double candidate[3];
            std::copy(origin, origin + 3, candidate);
            if (direction < 6)
            {
                candidate[direction / 2] += (direction % 2 ? 1 : -1) * radius * scale;
            }
            else if (direction < 9)
            {
                for (int i = 0; i < 3; i++)
                {
                    candidate[i] += (center[i] - origin[i]) * (direction == 6 ? .5 : direction == 7 ? .25
                                                                                                    : .125);
                }
            }
            else if (direction < 21)
            {
                const int plane = (direction - 9) / 4, signs = (direction - 9) % 4;
                const int a = plane == 2 ? 1 : 0, b = plane == 0 ? 1 : 2;
                candidate[a] += (signs & 1 ? 1 : -1) * radius * scale / std::sqrt(2.0);
                candidate[b] += (signs & 2 ? 1 : -1) * radius * scale / std::sqrt(2.0);
            }
            else
            {
                for (int i = 0; i < 3; i++)
                {
                    candidate[i] += ((direction - 21) & (1 << i) ? 1 : -1) * radius * scale / std::sqrt(3.0);
                }
            }
            // Clip the trial step to the exact positive-orientation kernel of
            // this star. Thin stars otherwise reject every fixed-size trial.
            double fraction = 1.0;
            for (int cell : star)
            {
                const auto &v = mesh.getTetra(cell);
                const double *start[4], *end[4];
                for (int i = 0; i < 4; i++)
                {
                    start[i] = v.v[i] == id ? origin : mesh.getPoint(v.v[i]);
                    end[i] = v.v[i] == id ? candidate : mesh.getPoint(v.v[i]);
                }
                const double a = RGeometricPredicates::orient3d(start[0], start[1], start[2], start[3]);
                const double b = RGeometricPredicates::orient3d(end[0], end[1], end[2], end[3]);
                if (b < a)
                {
                    fraction = std::min(fraction, .8 * a / (a - b));
                }
            }
            for (int i = 0; i < 3; i++)
            {
                candidate[i] = origin[i] + fraction * (candidate[i] - origin[i]);
            }
            const auto value = objective(candidate);
            if (improves(value, bestValue))
            {
                bestValue = value;
                std::copy(candidate, candidate + 3, best);
            }
        }
    }
    if (!improves(bestValue, std::array<double, 2>{std::min(angleBefore, std::max(1.0, options.getMinDihedralAngle())), before}))
    {
        return false;
    }
    mesh.setPoint(id, best);
    changed = star;
    touch(changed);
    geometryRevision[size_t(id)] = revision;
    profile.acceptedChanges++;
    return true;
}

bool RTetraMeshOptimizer::centroidRelocate(int id)
{
    RTetraMeshProfile::Timer timer(profile.relocation);
    if (boundary.getVertexType(id) != RBoundaryRecovery::VertexVolume || boundary.isBoundaryVertex(id))
    {
        return false;
    }
    std::vector<int> star;
    findStar(id, star);
    if (star.empty())
    {
        return false;
    }
    std::set<int> neighbors;
    double before = 1, maxVolume = 0, sizeBefore = 0;
    for (int cell : star)
    {
        if (mesh.isGhost(cell) || mesh.getTetra(cell).region < 0)
        {
            return false;
        }
        const auto &v = mesh.getTetra(cell);
        for (int p : v.v)
        {
            if (p != id)
            {
                neighbors.insert(p);
            }
        }
        const auto q = element(mesh, boundary, v.v);
        before = std::min(before, objectiveScore(q));
        maxVolume = std::max(maxVolume, q.volume);
        if (options.getUseSizeFunction())
        {
            sizeBefore = std::max(sizeBefore, sizeViolation(boundary, v.v, q));
        }
    }

    struct Cell
    {
        std::array<const double *, 4> points;
        std::array<bool, 4> fixed;
        std::array<int, 4> vertices;
        int moving;
    };

    std::vector<Cell> cells;
    cells.reserve(star.size());
    for (int idT : star)
    {
        const auto &t = mesh.getTetra(idT);
        Cell cell;
        for (int i = 0; i < 4; i++)
        {
            cell.points[i] = mesh.getPoint(t.v[i]);
            cell.vertices[i] = t.v[i];
            if (t.v[i] == id)
            {
                cell.moving = i;
            }
            const int *f = RTetraMesh::faceVertices[i];
            cell.fixed[i] = boundary.findSubFace(t.v[f[0]], t.v[f[1]], t.v[f[2]]) >= 0;
        }
        cells.push_back(cell);
    }
    double old[3], target[3] = {};
    std::copy(mesh.getPoint(id), mesh.getPoint(id) + 3, old);
    for (int v : neighbors)
    {
        for (int i = 0; i < 3; i++)
        {
            target[i] += mesh.getPoint(v)[i] / double(neighbors.size());
        }
    }
    for (double step : {.5, .25, .125, .0625})
    {
        double p[3];
        for (int i = 0; i < 3; i++)
        {
            p[i] = old[i] + step * (target[i] - old[i]);
        }
        bool valid = true;
        double after = 1;
        for (const auto &cell : cells)
        {
            auto points = cell.points;
            points[size_t(cell.moving)] = p;
            if (!(RGeometricPredicates::orient3d(points[0], points[1], points[2], points[3]) > 0))
            {
                valid = false;
                break;
            }
            profile.qualityMeasurements++;
            const auto q = RTetraMeshQuality::measure(points[0], points[1], points[2], points[3], cell.fixed.data());
            if (!(q.volume > 0) || (options.getVolumeConstraint() > 0 && q.volume > std::max(maxVolume, options.getVolumeConstraint()) * (1 + 1e-10)) || (options.getUseSizeFunction() && sizeViolation(boundary, cell.vertices.data(), q) > std::max(1.0, sizeBefore) * (1 + 1e-10)))
            {
                valid = false;
                break;
            }
            after = std::min(after, objectiveScore(q));
            if (!(after > before * (1 + 1e-5) + 1e-12))
            {
                valid = false;
                break;
            }
        }
        if (valid && after > before * (1 + 1e-5) + 1e-12)
        {
            mesh.setPoint(id, p);
            changed = star;
            touch(changed);
            geometryRevision[size_t(id)] = revision;
            profile.acceptedChanges++;
            return true;
        }
    }
    return false;
}

bool RTetraMeshOptimizer::placePoint(int id)
{
    if (options.getUseSizeFunction())
    {
        rejected.placement++;
        return false;
    }
    if (insertedPoints >= 16 || (options.getMaxSteinerNodes() > 0 && boundary.getNSteinerVertices() >= options.getMaxSteinerNodes()))
    {
        return false;
    }
    const auto t = mesh.getTetra(id);
    if (t.region < 0 || mesh.isGhost(id))
    {
        return false;
    }
    double center[3] = {0, 0, 0}, incenter[3] = {0, 0, 0};
    double areas[4], sum = 0;
    for (int i = 0; i < 4; i++)
    {
        const int *f = RTetraMesh::faceVertices[i];
        const double *a = mesh.getPoint(t.v[f[0]]), *b = mesh.getPoint(t.v[f[1]]), *c = mesh.getPoint(t.v[f[2]]);
        double u[3], v[3];
        for (int k = 0; k < 3; k++)
        {
            u[k] = b[k] - a[k];
            v[k] = c[k] - a[k];
        }
        const double n[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
        areas[i] = std::hypot(n[0], n[1], n[2]);
        sum += areas[i];
        for (int k = 0; k < 3; k++)
        {
            center[k] += .25 * mesh.getPoint(t.v[i])[k];
        }
    }
    if (!(sum > 0))
    {
        return false;
    }
    double shortest = std::numeric_limits<double>::max(), longest = 0;
    for (int i = 0; i < 4; i++)
    {
        for (int k = 0; k < 3; k++)
        {
            incenter[k] += areas[i] / sum * mesh.getPoint(t.v[i])[k];
        }
        for (int j = i + 1; j < 4; j++)
        {
            const double *a = mesh.getPoint(t.v[i]), *b = mesh.getPoint(t.v[j]);
            const double length = std::hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
            shortest = std::min(shortest, length);
            longest = std::max(longest, length);
        }
    }
    const double minRadius = std::max(.5 * shortest, 2e-8 * longest);
    const auto source = element(mesh, boundary, t.v);
    std::vector<std::array<double, 3>> candidates;
    for (double blend : {1.0, .5, 0.0})
    {
        std::array<double, 3> p;
        for (int k = 0; k < 3; k++)
        {
            p[k] = blend * incenter[k] + (1 - blend) * center[k];
        }
        candidates.push_back(p);
    }
    for (double blend : {.25, .5, 1.0})
    {
        std::array<double, 3> p;
        for (int k = 0; k < 3; k++)
        {
            p[k] = center[k] + blend * (source.circumcenter[size_t(k)] - center[k]);
        }
        candidates.push_back(p);
    }
    // Test off-centers in a bounded neighborhood. The constrained cavity must
    // contain the source cell; a candidate across a wall cannot repair this patch.
    for (const auto &candidate : candidates)
    {
        const double *p = candidate.data();
        if (!std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2]) || std::hypot(p[0] - center[0], p[1] - center[1], p[2] - center[2]) > 2 * longest)
        {
            rejected.placement++;
            continue;
        }
        const int seed = mesh.locate(p, id);
        if (mesh.isGhost(seed) || mesh.getTetra(seed).region != t.region || mesh.findDuplicateVertex(seed, p) >= 0 || !boundary.prepareVolumeCavity(p, seed))
        {
            rejected.placement++;
            continue;
        }
        if (std::find(mesh.getCavityTetras().begin(), mesh.getCavityTetras().end(), id) == mesh.getCavityTetras().end())
        {
            rejected.placement++;
            continue;
        }
        if (boundary.findEncroachment(p).type != RBoundaryRecovery::EncroachmentNone)
        {
            rejected.boundary++;
            continue;
        }
        if (!boundary.insertionAllowed(p, RBoundaryRecovery::VertexVolume, minRadius))
        {
            rejected.spacing++;
            continue;
        }
        double before = 1, after = 1, maxVolume = 0, angleBefore = 180, angleAfter = 180, ratioBefore = 0, ratioAfter = 0;
        bool valid = true;
        size_t subBefore = 0, subAfter = 0;
        for (int cell : mesh.getCavityTetras())
        {
            if (mesh.isGhost(cell) || mesh.getTetra(cell).region != t.region)
            {
                valid = false;
                break;
            }
            const auto &v = mesh.getTetra(cell);
            const auto q = element(mesh, boundary, v.v);
            before = std::min(before, objectiveScore(q));
            maxVolume = std::max(maxVolume, q.volume);
            subBefore += q.minFreeAngle < 1;
            angleBefore = std::min(angleBefore, q.minFreeAngle);
            ratioBefore = std::max(ratioBefore, q.radiusEdge);
        }
        for (const auto &face : mesh.getCavityFaces())
        {
            const auto &v = mesh.getTetra(face.tetra);
            const int *f = RTetraMesh::faceVertices[face.face];
            if (v.v[f[0]] < 0 || v.v[f[1]] < 0 || v.v[f[2]] < 0)
            {
                valid = false;
                break;
            }
            const auto q = RTetraMeshQuality::measure(mesh.getPoint(v.v[f[0]]), mesh.getPoint(v.v[f[1]]), mesh.getPoint(v.v[f[2]]), p);
            if (!(q.volume > 0) || (options.getVolumeConstraint() > 0 && q.volume > std::max(maxVolume, options.getVolumeConstraint()) * (1 + 1e-10)))
            {
                valid = false;
                break;
            }
            after = std::min(after, objectiveScore(q));
            subAfter += q.minFreeAngle < 1;
            angleAfter = std::min(angleAfter, q.minAngle);
            ratioAfter = std::max(ratioAfter, q.radiusEdge);
        }
        if (angleAfter + 1e-8 < std::min(angleBefore, 1.0) || ratioAfter > std::max(ratioBefore, options.getRadiusEdgeRatio()) * (1 + 1e-8))
        {
            valid = false;
        }
        if (!valid || subAfter > subBefore || !(after > before * 1.05 + 1e-12))
        {
            rejected.placement++;
            continue;
        }
        boundary.insertCavityPoint(p, RBoundaryRecovery::VertexVolume, minRadius);
        changed = mesh.getCreatedTetras();
        touch(changed);
        boundaryEpoch++;
        insertedPoints++;
        profile.acceptedChanges++;
        return true;
    }
    return false;
}

bool RTetraMeshOptimizer::remeshPatch(int id)
{
    if (remeshAttempts++ >= 128)
    {
        return false;
    }
    const int region = mesh.getTetra(id).region;
    std::set<int> cells = {id};
    // Grow only across unconstrained faces until the cavity is convex. The
    // complete local triangulation is then committed atomically, never as an
    // intermediate sequence of inverted or boundary-breaking flips.
    for (int step = 0; step < 24; step++)
    {
        std::set<int> points;
        for (int t : cells)
        {
            for (int v : mesh.getTetra(t).v)
            {
                points.insert(v);
            }
        }
        if (points.size() > 20 || cells.size() > 24)
        {
            rejected.topology++;
            return false;
        }
        int expand = -1;
        for (int t : cells)
        {
            const auto &tet = mesh.getTetra(t);
            for (int face = 0; face < 4; face++)
            {
                if (!cells.count(tet.n[face]))
                {
                    const int *f = RTetraMesh::faceVertices[face];
                    bool convex = true;
                    for (int v : points)
                    {
                        if (RGeometricPredicates::orient3d(mesh.getPoint(tet.v[f[0]]), mesh.getPoint(tet.v[f[1]]), mesh.getPoint(tet.v[f[2]]), mesh.getPoint(v)) < 0)
                        {
                            convex = false;
                            break;
                        }
                    }
                    if (!convex)
                    {
                        if (boundary.findSubFace(tet.v[f[0]], tet.v[f[1]], tet.v[f[2]]) >= 0 || mesh.isGhost(tet.n[face]) || mesh.getTetra(tet.n[face]).region != region)
                        {
                            rejected.boundary++;
                            return false;
                        }
                        expand = tet.n[face];
                        break;
                    }
                }
            }
            if (expand >= 0)
            {
                break;
            }
        }
        if (expand >= 0)
        {
            cells.insert(expand);
            continue;
        }
        // A single tetrahedron cannot be improved by retriangulating its vertices.
        if (cells.size() == 1)
        {
            const auto &t = mesh.getTetra(id);
            for (int face = 0; face < 4; face++)
            {
                const int *f = RTetraMesh::faceVertices[face];
                const int n = t.n[face];
                if (!mesh.isGhost(n) && mesh.getTetra(n).region == region && boundary.findSubFace(t.v[f[0]], t.v[f[1]], t.v[f[2]]) < 0)
                {
                    cells.insert(n);
                }
            }
            if (cells.size() == 1)
            {
                return false;
            }
            continue;
        }
        RTetraMesh local;
        std::vector<int> map(points.begin(), points.end());
        for (int v : map)
        {
            const double *p = mesh.getPoint(v);
            local.addPoint(p[0], p[1], p[2]);
        }
        std::vector<int> duplicates;
        try
        {
            RDelaunayTetrahedralizer::tetrahedralize(local, duplicates);
        }
        catch (const RError &)
        {
            rejected.topology++;
            return false;
        }
        if (std::any_of(duplicates.begin(), duplicates.end(), [](int v)
                        {
                            return v >= 0;
                        }))
        {
            rejected.topology++;
            return false;
        }
        std::vector<std::array<int, 4>> proposed;
        for (int t = 0; t < local.getNTetraSlots(); t++)
        {
            if (local.isAlive(t) && !local.isGhost(t))
            {
                const auto &v = local.getTetra(t);
                proposed.push_back({map[size_t(v.v[0])], map[size_t(v.v[1])], map[size_t(v.v[2])], map[size_t(v.v[3])]});
            }
        }
        return replace(std::vector<int>(cells.begin(), cells.end()), proposed, false);
    }
    rejected.topology++;
    return false;
}

bool RTetraMeshOptimizer::compoundFlip(int id)
{
    RTetraMeshProfile::Timer timer(profile.compound);
    if (compoundAttempts++ >= 32)
    {
        return false;
    }

    struct Proposal
    {
        std::vector<int> removed;
        std::vector<std::array<int, 4>> added;
    };

    auto orient = [&](RTetraMesh &trial, Proposal &p)
    {
        for (auto &v : p.added)
        {
            const double sign = RGeometricPredicates::orient3d(trial.getPoint(v[0]), trial.getPoint(v[1]), trial.getPoint(v[2]), trial.getPoint(v[3]));
            if (sign == 0)
            {
                return false;
            }
            if (sign < 0)
            {
                std::swap(v[0], v[1]);
            }
        }
        return true;
    };
    auto faceProposal = [&](RTetraMesh &trial, const std::array<int, 3> &key, Proposal &p)
    {
        if (boundary.findSubFace(key[0], key[1], key[2]) >= 0)
        {
            return false;
        }
        RTetraMesh::Face face;
        if (!trial.findFace(key[0], key[1], key[2], face))
        {
            return false;
        }
        const auto &t = trial.getTetra(face.tetra);
        const int n = t.n[face.face];
        if (trial.isGhost(face.tetra) || trial.isGhost(n) || t.region != trial.getTetra(n).region)
        {
            return false;
        }
        const int d = t.v[face.face];
        int e = -1;
        for (int v : trial.getTetra(n).v)
        {
            if (v != key[0] && v != key[1] && v != key[2])
            {
                e = v;
            }
        }
        if (e < 0 || trial.findEdge(d, e) != RTetraMesh::noTetra)
        {
            return false;
        }
        p.removed = {face.tetra, n};
        p.added = {{d, e, key[0], key[1]}, {d, e, key[1], key[2]}, {d, e, key[2], key[0]}};
        return orient(trial, p);
    };
    const auto seed = mesh.getTetra(id);
    for (int face = 0; face < 4; face++)
    {
        const int *f = RTetraMesh::faceVertices[face];
        const std::array<int, 3> key = {seed.v[f[0]], seed.v[f[1]], seed.v[f[2]]};
        RTetraMesh::ReplacementTrial firstTrial(mesh);
        RTetraMesh &first = mesh;
        Proposal p;
        if (!faceProposal(first, key, p) || !firstTrial.replace(p.removed, p.added))
        {
            continue;
        }
        profile.trialReplacements++;
        const auto firstAdded = first.getCreatedTetras();
        std::set<std::array<int, 3>> faces;
        std::set<std::array<int, 2>> edges;
        for (int t : firstAdded)
        {
            const auto &v = first.getTetra(t);
            for (int i = 0; i < 4; i++)
            {
                const int *side = RTetraMesh::faceVertices[i];
                auto k = std::array<int, 3>{v.v[side[0]], v.v[side[1]], v.v[side[2]]};
                std::sort(k.begin(), k.end());
                faces.insert(k);
                for (int j = i + 1; j < 4; j++)
                {
                    edges.insert({std::min(v.v[i], v.v[j]), std::max(v.v[i], v.v[j])});
                }
            }
        }
        std::vector<Proposal> seconds;
        for (const auto &k : faces)
        {
            Proposal next;
            if (faceProposal(first, k, next))
            {
                seconds.push_back(next);
            }
        }
        for (const auto &edge : edges)
        {
            if (boundary.findSubSegment(edge[0], edge[1]) >= 0)
            {
                continue;
            }
            std::vector<int> star;
            first.findStar(edge[0], star);
            Proposal next;
            std::set<int> ring;
            bool valid = true;
            for (int t : star)
            {
                const auto &v = first.getTetra(t);
                if (std::find(v.v, v.v + 4, edge[1]) == v.v + 4)
                {
                    continue;
                }
                if (first.isGhost(t))
                {
                    valid = false;
                    break;
                }
                next.removed.push_back(t);
                for (int point : v.v)
                {
                    if (point != edge[0] && point != edge[1])
                    {
                        ring.insert(point);
                    }
                }
            }
            if (!valid || next.removed.size() != 3 || ring.size() != 3)
            {
                continue;
            }
            for (int v : ring)
            {
                if (boundary.findSubFace(edge[0], edge[1], v) >= 0)
                {
                    valid = false;
                }
            }
            if (!valid)
            {
                continue;
            }
            auto it = ring.begin();
            const int x = *it++, y = *it++, z = *it;
            next.added = {{x, y, z, edge[0]}, {x, z, y, edge[1]}};
            if (orient(first, next))
            {
                seconds.push_back(next);
            }
        }
        std::vector<Proposal> complete;
        for (const auto &next : seconds)
        {
            RTetraMesh::ReplacementTrial secondTrial(mesh);
            RTetraMesh &second = mesh;
            if (!secondTrial.replace(next.removed, next.added))
            {
                continue;
            }
            profile.trialReplacements++;
            profile.journalBytes += secondTrial.bytesSaved();
            std::set<int> removed(p.removed.begin(), p.removed.end()), added;
            for (int old : next.removed)
            {
                if (std::find(firstAdded.begin(), firstAdded.end(), old) == firstAdded.end())
                {
                    removed.insert(old);
                }
            }
            for (int cell : firstAdded)
            {
                if (second.isAlive(cell))
                {
                    added.insert(cell);
                }
            }
            added.insert(second.getCreatedTetras().begin(), second.getCreatedTetras().end());
            std::vector<std::array<int, 4>> proposed;
            for (int cell : added)
            {
                const auto &v = second.getTetra(cell);
                proposed.push_back({v.v[0], v.v[1], v.v[2], v.v[3]});
            }
            complete.push_back({std::vector<int>(removed.begin(), removed.end()), std::move(proposed)});
        }
        profile.journalBytes += firstTrial.bytesSaved();
        firstTrial.rollback();
        // Restore the original topology before evaluating each complete patch.
        // Candidate and slot order match the former independent trial copies.
        for (auto &candidate : complete)
        {
            if (replace(candidate.removed, candidate.added, false))
            {
                return true;
            }
        }
    }
    return false;
}

bool RTetraMeshOptimizer::boundaryPatch(int id)
{
    RTetraMeshProfile::Timer timer(profile.boundaryTrials);
    if (boundaryAttempts >= 24 || options.getUseSizeFunction())
    {
        return false;
    }
    const auto seed = mesh.getTetra(id);
    std::vector<int> faces;
    for (int i = 0; i < 4; i++)
    {
        const int *f = RTetraMesh::faceVertices[i];
        const int face = boundary.findSubFace(seed.v[f[0]], seed.v[f[1]], seed.v[f[2]]);
        if (face >= 0)
        {
            faces.push_back(face);
        }
    }
    if (faces.empty())
    {
        return false;
    }

    struct Metrics
    {
        size_t subDegree = 0, badAngle = 0, nearFlat = 0;
        double minFree = 180, minAngle = 180, volume = 0, area = 0, maxVolume = 0, maxRadius = 0, radius95 = 0;
    };

    auto evaluate = [&](const RTetraMesh &m, const RBoundaryRecovery &b)
    {
        Metrics result;
        std::vector<double> radii;
        for (int t = 0; t < m.getNTetraSlots(); t++)
        {
            if (m.isAlive(t) && !m.isGhost(t) && m.getTetra(t).region >= 0)
            {
                const auto q = element(m, b, m.getTetra(t).v);
                result.subDegree += q.minFreeAngle < 1;
                result.badAngle += q.minFreeAngle < std::max(1.0, options.getMinDihedralAngle());
                result.minFree = std::min(result.minFree, q.minFreeAngle);
                result.minAngle = std::min(result.minAngle, q.minAngle);
                result.volume += q.volume;
                result.maxVolume = std::max(result.maxVolume, q.volume);
                result.maxRadius = std::max(result.maxRadius, q.radiusEdge);
                radii.push_back(q.radiusEdge);
                if (q.scaledVolume < 1e-12 || !(q.volume > 0))
                {
                    result.nearFlat++;
                }
            }
        }
        for (const auto &face : b.getSubFaces())
        {
            if (face.alive)
            {
                const double *a = m.getPoint(face.v[0]), *p = m.getPoint(face.v[1]), *q = m.getPoint(face.v[2]);
                const double u[3] = {p[0] - a[0], p[1] - a[1], p[2] - a[2]}, v[3] = {q[0] - a[0], q[1] - a[1], q[2] - a[2]};
                result.area += .5 * std::hypot(u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]);
            }
        }
        std::sort(radii.begin(), radii.end());
        if (!radii.empty())
        {
            result.radius95 = radii[size_t(.95 * double(radii.size() - 1))];
        }
        return result;
    };
    const auto before = evaluate(mesh, boundary);
    std::set<int> regions;
    for (int t = 0; t < mesh.getNTetraSlots(); t++)
    {
        if (mesh.isAlive(t) && mesh.getTetra(t).region >= 0)
        {
            regions.insert(mesh.getTetra(t).region);
        }
    }
    double shortest = std::numeric_limits<double>::max();
    for (int i = 0; i < 4; i++)
    {
        for (int j = i + 1; j < 4; j++)
        {
            const double *a = mesh.getPoint(seed.v[i]), *b = mesh.getPoint(seed.v[j]);
            shortest = std::min(shortest, std::hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]));
        }
    }
    for (int face : faces)
    {
        for (bool centroid : {false, true})
        {
            if (boundaryAttempts++ >= 24)
            {
                return false;
            }
            profile.trialCopies++;
            profile.trialCopyBytes += mesh.storageBytes();
            RTetraMesh trial = mesh;
            RBoundaryRecovery trialBoundary(trial, boundary);
            const auto count = boundary.getNSteinerVertices();
            const auto userLimit = options.getMaxSteinerNodes() > 0 ? uint64_t(options.getMaxSteinerNodes()) : count + 8;
            trialBoundary.setMaxSteinerVertices(std::min(count + 8, userLimit));
            try
            {
                if (!trialBoundary.splitSubFaceAt(face, .25 * shortest, centroid))
                {
                    rejected.spacing++;
                    continue;
                }
                trialBoundary.recover();
                if (trialBoundary.classifyRegions() != int(regions.size()))
                {
                    rejected.region++;
                    continue;
                }
            }
            catch (const RError &)
            {
                rejected.placement++;
                continue;
            }
            const auto after = evaluate(trial, trialBoundary);
            const bool improves = after.subDegree < before.subDegree || (after.subDegree == before.subDegree && after.badAngle < before.badAngle);
            if (!improves || after.minFree + 1e-8 < before.minFree || after.minAngle + 1e-8 < before.minAngle || after.nearFlat > before.nearFlat || after.maxRadius > before.maxRadius * (1 + 1e-8) || after.radius95 > before.radius95 * (1 + 1e-8) || (options.getVolumeConstraint() > 0 && after.maxVolume > std::max(before.maxVolume, options.getVolumeConstraint()) * (1 + 1e-10)) || std::fabs(after.volume - before.volume) > 1e-10 * std::max(before.volume, 1e-30) || std::fabs(after.area - before.area) > 1e-10 * std::max(before.area, 1e-30))
            {
                rejected.quality++;
                continue;
            }
            // Trial recovery only subdivides existing facets/segments; no original
            // coordinate is moved. Reject any change to an existing coordinate anyway.
            bool fixed = true;
            for (int v = 0; v < mesh.getNPoints(); v++)
            {
                for (int i = 0; i < 3; i++)
                {
                    fixed = fixed && trial.getPoint(v)[i] == mesh.getPoint(v)[i];
                }
            }
            if (!fixed)
            {
                rejected.boundary++;
                continue;
            }
            changed.clear();
            for (int t = 0; t < trial.getNTetraSlots(); t++)
            {
                if (trial.isAlive(t) && !trial.isGhost(t))
                {
                    const auto &v = trial.getTetra(t);
                    if (t >= mesh.getNTetraSlots() || !mesh.isAlive(t) || v.region != mesh.getTetra(t).region || !std::equal(v.v, v.v + 4, mesh.getTetra(t).v) || !std::equal(v.n, v.n + 4, mesh.getTetra(t).n))
                    {
                        changed.push_back(t);
                    }
                }
            }
            boundaryPoints += size_t(trialBoundary.getNSteinerVertices() - count);
            trialBoundary.setMaxSteinerVertices(boundary.maxSteinerVertices);
            mesh = std::move(trial);
            boundary.copyState(trialBoundary);
            touch(changed);
            boundaryEpoch++;
            profile.acceptedChanges++;
            return true;
        }
    }
    return false;
}

unsigned int RTetraMeshOptimizer::orderedRepair(bool moveVertices)
{
    RTetraMeshProfile::Timer timer(profile.orderedRepair);
    RTetraMeshProgress progress(moveVertices ? "Optimizing mesh: preparing repair sweeps" : "Repairing initial mesh");
    orderedRepairActive = true;
    unsigned int changes = 0;
    // Vertices of each cell when it was last scheduled, and the list generation it went to.
    // A cell met again unchanged would only add keys which its list already holds.
    std::vector<std::array<int, 4>> scheduledCells;
    std::vector<unsigned int> scheduledMarks;
    unsigned int scheduleMark = 0;
    for (int stage = 0; stage < 2; stage++)
    {
        progress.initialize(stage == 0 ? "Repairing mesh: preparing Delaunay sweeps" : "Optimizing mesh: preparing quality sweeps");
        // Work for the next sweep is appended as it is found and put in order when the
        // sweep starts, which visits the same keys in the same order as ordered sets would.
        std::vector<std::array<int, 3>> faces;
        std::vector<std::array<int, 2>> edges;
        std::vector<int> vertices;
        size_t faceLimit = 0, edgeLimit = 0, vertexLimit = 0;
        auto schedule = [&](int id)
        {
            if (!mesh.isAlive(id) || mesh.isGhost(id))
            {
                return;
            }
            const auto &t = mesh.getTetra(id);
            if (size_t(id) >= scheduledMarks.size())
            {
                scheduledMarks.resize(size_t(mesh.getNTetraSlots()), 0);
                scheduledCells.resize(scheduledMarks.size());
            }
            if (scheduledMarks[size_t(id)] == scheduleMark && std::equal(t.v, t.v + 4, scheduledCells[size_t(id)].begin()))
            {
                return;
            }
            scheduledMarks[size_t(id)] = scheduleMark;
            std::copy(t.v, t.v + 4, scheduledCells[size_t(id)].begin());
            if (moveVertices && stage == 1)
            {
                for (int v : t.v)
                {
                    if (boundary.getVertexType(v) == RBoundaryRecovery::VertexVolume && !boundary.isBoundaryVertex(v))
                    {
                        vertices.push_back(v);
                    }
                }
            }
            for (int i = 0; i < 4; i++)
            {
                const int *f = RTetraMesh::faceVertices[i];
                auto key = std::array<int, 3>{t.v[f[0]], t.v[f[1]], t.v[f[2]]};
                std::sort(key.begin(), key.end());
                faces.push_back(key);
                for (int j = i + 1; j < 4; j++)
                {
                    edges.push_back({std::min(t.v[i], t.v[j]), std::max(t.v[i], t.v[j])});
                }
            }
        };
        // Repeated keys are dropped whenever a list doubles, which bounds its memory.
        auto compact = [](auto &keys, size_t &limit)
        {
            if (keys.size() < limit)
            {
                return;
            }
            makeUnique(keys);
            limit = std::max(size_t(65536), 2 * keys.size());
        };
        faces.reserve(4 * size_t(mesh.getNTetraSlots()));
        edges.reserve(6 * size_t(mesh.getNTetraSlots()));
        ++scheduleMark;
        for (int id = 0; id < mesh.getNTetraSlots(); id++)
        {
            schedule(id);
        }
        for (int wave = 0; wave < (stage == 0 ? 5 : 8) && (!faces.empty() || !edges.empty() || !vertices.empty()); wave++)
        {
            makeUnique(faces);
            makeUnique(edges);
            makeUnique(vertices);
            const auto currentFaces = std::move(faces);
            const auto currentEdges = std::move(edges);
            const auto currentVertices = std::move(vertices);
            faces.clear();
            edges.clear();
            vertices.clear();
            faceLimit = edgeLimit = vertexLimit = 65536;
            ++scheduleMark;
            progress.initialize(QString("%1: sweep %2 of %3")
                                    .arg(stage == 0 ? "Repairing mesh Delaunayness" : "Optimizing mesh quality")
                                    .arg(wave + 1)
                                    .arg(stage == 0 ? 5 : 8),
                                false);
            const uint64_t total = currentFaces.size() + currentEdges.size() + currentVertices.size();
            uint64_t processed = 0;
            auto reschedule = [&]()
            {
                std::set<int> points;
                for (int id : changed)
                {
                    if (mesh.isAlive(id))
                    {
                        for (int v : mesh.getTetra(id).v)
                        {
                            points.insert(v);
                        }
                    }
                }
                std::vector<int> star, affected;
                for (int v : points)
                {
                    findStar(v, star);
                    affected.insert(affected.end(), star.begin(), star.end());
                }
                std::sort(affected.begin(), affected.end());
                affected.erase(std::unique(affected.begin(), affected.end()), affected.end());
                // Each cell contributes the same set keys once, even when it
                // belongs to several changed vertex stars.
                for (int id : affected)
                {
                    schedule(id);
                }
                compact(faces, faceLimit);
                compact(edges, edgeLimit);
                compact(vertices, vertexLimit);
            };
            for (const auto &f : currentFaces)
            {
                orderedVisits++;
                if (flipFace(f[0], f[1], f[2], stage == 0))
                {
                    changes++;
                    reschedule();
                }
                progress.print(++processed, total);
            }
            for (const auto &e : currentEdges)
            {
                orderedVisits++;
                if (flipEdge(e[0], e[1], stage == 0))
                {
                    changes++;
                    reschedule();
                }
                progress.print(++processed, total);
            }
            for (int v : currentVertices)
            {
                orderedVisits++;
                if (centroidRelocate(v))
                {
                    changes++;
                    reschedule();
                }
                progress.print(++processed, total);
            }
        }
        if (!faces.empty() || !edges.empty() || !vertices.empty())
        {
            budgetStops++;
        }
    }
    orderedRepairActive = false;
    return changes;
}

unsigned int RTetraMeshOptimizer::optimize(bool smoothVertices)
{
    edgeRevision.clear();
    faceRevision.clear();
    // The caller may have refined or recovered the mesh between calls.
    qualityCache.clear();
    starCache.clear();
    starCacheStamp.clear();
    vertexRevision.assign(size_t(mesh.getNPoints()), ++revision);
    geometryRevision = vertexRevision;
    boundaryEpoch++;
    if (!smoothVertices)
    {
        return orderedRepair();
    }

    struct Candidate
    {
        double priority;
        int id;
        std::array<int, 4> vertices;

        bool operator<(const Candidate &b) const
        {
            return priority == b.priority ? vertices > b.vertices : priority < b.priority;
        }
    };

    unsigned int changes = orderedRepair(true);
    RTetraMeshProgress progress("Optimizing mesh: preparing bounded patch repairs");
    remeshAttempts = compoundAttempts = boundaryAttempts = 0;
    // Ordered repair already restored Delaunayness and quality. The finishing
    // queue must not run another Delaunay stage which can undo sliver repairs.
    for (int stage = 1; stage < 2; stage++)
    {
        std::priority_queue<Candidate> queue;
        std::unordered_map<std::array<int, 4>, int, KeyHash> pending;
        std::unordered_map<std::array<int, 4>, unsigned int, KeyHash> visits;
        std::vector<unsigned int> vertexVisits(size_t(mesh.getNPoints()), 0);
        std::vector<size_t> vertexAttemptStamp(size_t(mesh.getNPoints()), 0);
        auto enqueue = [&](int id)
        {
            if (id < 0 || !mesh.isAlive(id) || mesh.isGhost(id) || (stage == 1 && mesh.getTetra(id).region < 0))
            {
                return;
            }
            const auto &t = mesh.getTetra(id);
            std::array<int, 4> key = {t.v[0], t.v[1], t.v[2], t.v[3]};
            std::sort(key.begin(), key.end());
            if (visits[key] >= 8)
            {
                return;
            }
            auto queued = pending.find(key);
            if (queued != pending.end() && queued->second == id)
            {
                return;
            }
            pending[key] = id;
            const auto q = element(mesh, boundary, t.v);
            if (stage == 1 && q.radiusEdge <= options.getRadiusEdgeRatio() && q.minFreeAngle >= std::max(10.0, options.getMinDihedralAngle()) && q.scaledVolume >= 1e-12)
            {
                pending.erase(key);
                return;
            }
            queue.push({q.minFreeAngle < 1 ? 100 + (1 - q.minFreeAngle) : -objectiveScore(q), id, key});
        };
        for (int id = 0; id < mesh.getNTetraSlots(); id++)
        {
            enqueue(id);
        }
        const size_t budget = queue.size() * (stage == 0 ? 4 : 8) + 64;
        size_t processed = 0;
        while (!queue.empty() && processed < budget)
        {
            progress.pulse([&]()
                           {
                               return QString("Optimizing mesh: %1 patches visited, %2 queued")
                                   .arg(qulonglong(processed))
                                   .arg(qulonglong(queue.size()));
                           });
            const auto entry = queue.top();
            queue.pop();
            auto queued = pending.find(entry.vertices);
            if (queued == pending.end() || queued->second != entry.id)
            {
                continue;
            }
            pending.erase(queued);
            if (!mesh.isAlive(entry.id) || mesh.isGhost(entry.id))
            {
                continue;
            }
            const auto t = mesh.getTetra(entry.id);
            auto key = std::array<int, 4>{t.v[0], t.v[1], t.v[2], t.v[3]};
            std::sort(key.begin(), key.end());
            if (key != entry.vertices || visits[key] >= 8)
            {
                continue;
            }
            visits[key]++;
            processed++;
            queueVisits++;
            const auto q = element(mesh, boundary, t.v);
            if (stage == 1 && q.radiusEdge <= options.getRadiusEdgeRatio() && q.minFreeAngle >= std::max(10.0, options.getMinDihedralAngle()) && q.scaledVolume >= 1e-12)
            {
                continue;
            }
            const auto previous = rejected;
            changed.clear();
            bool repaired = false;
            for (int i = 0; i < 4 && !repaired; i++)
            {
                const int *f = RTetraMesh::faceVertices[i];
                repaired = flipFace(t.v[f[0]], t.v[f[1]], t.v[f[2]], stage == 0);
            }
            for (int i = 0; i < 6 && !repaired; i++)
            {
                repaired = flipEdge(t.v[RTetraMesh::edgeVertices[i][0]], t.v[RTetraMesh::edgeVertices[i][1]], stage == 0);
            }
            if (!repaired && stage == 1 && q.scaledVolume < 1e-12)
            {
                repaired = compoundFlip(entry.id);
                if (!repaired)
                {
                    repaired = remeshPatch(entry.id);
                }
            }
            if (!repaired && stage == 1 && smoothVertices)
            {
                vertexVisits.resize(size_t(mesh.getNPoints()), 0);
                vertexAttemptStamp.resize(size_t(mesh.getNPoints()), 0);
                for (int v : t.v)
                {
                    if (vertexVisits[size_t(v)] >= 8 || boundary.getVertexType(v) != RBoundaryRecovery::VertexVolume || boundary.isBoundaryVertex(v))
                    {
                        continue;
                    }
                    std::vector<int> star;
                    findStar(v, star);
                    size_t stamp = 1;
                    for (int cell : star)
                    {
                        for (int point : mesh.getTetra(cell).v)
                        {
                            if (point >= 0)
                            {
                                stamp = std::max(stamp, vertexRevision[size_t(point)]);
                            }
                        }
                    }
                    if (vertexAttemptStamp[size_t(v)] == stamp)
                    {
                        profile.queueSkips++;
                        continue;
                    }
                    vertexVisits[size_t(v)]++;
                    const bool moved = relocate(v);
                    vertexAttemptStamp[size_t(v)] = moved ? revision : stamp;
                    if (moved)
                    {
                        repaired = true;
                        break;
                    }
                }
                if (!repaired && (q.minFreeAngle < 1.0 || q.radiusEdge > 10.0))
                {
                    repaired = compoundFlip(entry.id);
                }
                if (!repaired && (q.minFreeAngle < 1.0 || q.radiusEdge > 10.0))
                {
                    repaired = remeshPatch(entry.id);
                }
                if (!repaired && q.minFreeAngle < 1.0)
                {
                    repaired = boundaryPatch(entry.id);
                }
                if (!repaired && (q.minFreeAngle < 5.0 || q.radiusEdge > std::max(3.0, 2.0 * options.getRadiusEdgeRatio())))
                {
                    repaired = placePoint(entry.id);
                }
            }
            auto &record = attempts[key];
            record.vertices = key;
            record.visits++;
            record.accepted += repaired;
            record.boundary += rejected.boundary - previous.boundary;
            record.region += rejected.region - previous.region;
            record.topology += rejected.topology - previous.topology;
            record.quality += rejected.quality - previous.quality;
            record.spacing += rejected.spacing - previous.spacing;
            record.placement += rejected.placement - previous.placement;
            if (!repaired)
            {
                continue;
            }
            changes++;
            // Edge-link legality can change beyond the face neighbors. Revisit
            // the complete vertex neighborhood, with canonical-cell deduplication.
            std::set<int> affected;
            for (int id : changed)
            {
                if (mesh.isAlive(id))
                {
                    for (int v : mesh.getTetra(id).v)
                    {
                        affected.insert(v);
                    }
                }
            }
            std::vector<int> star;
            for (int v : affected)
            {
                findStar(v, star);
                for (int id : star)
                {
                    enqueue(id);
                }
            }
        }
        if (!queue.empty())
        {
            budgetStops++;
        }
    }
    return changes;
}

std::vector<RTetraMeshQuality::PatchAttempts> RTetraMeshOptimizer::getPatchAttempts() const
{
    std::vector<RTetraMeshQuality::PatchAttempts> result;
    for (const auto &record : attempts)
    {
        result.push_back(record.second);
    }
    return result;
}

unsigned int RTetraMeshOptimizer::repairCreated(bool preferSimpleFlips)
{
    RTetraMeshProfile::Timer timer(profile.localRepair);
    simpleRepair = preferSimpleFlips;
    const auto seed = boundary.getCreatedTetras();
    boundaryEpoch++;
    touch(seed);
    std::vector<std::array<int, 3>> faces;
    std::vector<std::array<int, 2>> edges;
    faces.reserve(4 * seed.size());
    edges.reserve(6 * seed.size());
    for (int id : seed)
    {
        if (mesh.isAlive(id) && !mesh.isGhost(id))
        {
            const auto t = mesh.getTetra(id);
            for (int i = 0; i < 4; i++)
            {
                const int *f = RTetraMesh::faceVertices[i];
                std::array<int, 3> key = {t.v[f[0]], t.v[f[1]], t.v[f[2]]};
                std::sort(key.begin(), key.end());
                faces.push_back(key);
                for (int j = i + 1; j < 4; j++)
                {
                    edges.push_back({std::min(t.v[i], t.v[j]), std::max(t.v[i], t.v[j])});
                }
            }
        }
    }
    makeUnique(faces);
    makeUnique(edges);
    unsigned int count = 0;
    for (const auto &f : faces)
    {
        count += flipFace(f[0], f[1], f[2], true);
    }
    for (const auto &e : edges)
    {
        count += flipEdge(e[0], e[1], true);
    }
    simpleRepair = false;
    return count;
}
