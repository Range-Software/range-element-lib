#ifndef REL_TETRA_MESH_PROFILE_H
#define REL_TETRA_MESH_PROFILE_H

#include <chrono>
#include <cstdint>

//! Per-generation wall-clock phases and work counters; no global mutable state.
//! Optimizer timings are nested inside the generation phases, not additive to them.
struct RTetraMeshProfile
{
    //! Seconds; total excludes optional report serialization and final logs.
    double total = 0, tetrahedralization = 0, boundary = 0, initialRepair = 0, refinement = 0, optimization = 0, diagnostics = 0;
    double orderedRepair = 0, localRepair = 0, relocation = 0, compound = 0, boundaryTrials = 0;
    //! Optimizer work only. Copy bytes count mesh payload, excluding boundary metadata
    //! and allocator overhead; journal bytes count locally saved records.
    uint64_t qualityMeasurements = 0, qualityCacheHits = 0, starQueries = 0, starCacheHits = 0;
    uint64_t trialCopies = 0, trialCopyBytes = 0, journalBytes = 0, trialReplacements = 0, queueSkips = 0, acceptedChanges = 0;
    //! Delaunay flips rejected by an empty-sphere test before any proposal was built.
    uint64_t delaunaySkips = 0;

    void addWork(const RTetraMeshProfile &p)
    {
        orderedRepair += p.orderedRepair;
        localRepair += p.localRepair;
        relocation += p.relocation;
        compound += p.compound;
        boundaryTrials += p.boundaryTrials;
        qualityMeasurements += p.qualityMeasurements;
        qualityCacheHits += p.qualityCacheHits;
        starQueries += p.starQueries;
        starCacheHits += p.starCacheHits;
        trialCopies += p.trialCopies;
        trialCopyBytes += p.trialCopyBytes;
        journalBytes += p.journalBytes;
        trialReplacements += p.trialReplacements;
        queueSkips += p.queueSkips;
        acceptedChanges += p.acceptedChanges;
        delaunaySkips += p.delaunaySkips;
    }

    struct Timer
    {
        double &value;
        std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

        explicit Timer(double &value)
            : value(value)
        {
        }

        bool active = true;

        void finish()
        {
            if (active)
            {
                value += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
                active = false;
            }
        }

        ~Timer()
        {
            finish();
        }

        Timer(const Timer &) = delete;
    };
};
#endif
