#ifndef REL_TETRA_MESH_QUALITY_H
#define REL_TETRA_MESH_QUALITY_H

#include <cstddef>
#include <array>
#include <string>
#include <vector>

#include "rel_tetra_mesh_data.h"
#include "rel_tetra_mesh_options.h"
#include "rel_tetra_mesh_profile.h"

//! Scale independent element quality and requested-constraint diagnostics.
class RTetraMeshQuality
{
    public:

        struct Element
        {
            double volume = 0.0;
            double scaledVolume = 0.0; //!< Six times volume divided by longest edge cubed.
            double meanRatio = 0.0;    //!< One for a regular tetrahedron, zero for a flat one.
            double radiusEdge = 0.0;
            double circumradius = 0.0;
            std::array<double, 3> circumcenter = {};
            double minAngle = 0.0;
            double minFreeAngle = 0.0;
            double maxAngle = 180.0;
        };

        struct Summary
        {
            size_t elements = 0, invalid = 0, nearFlat = 0;
            size_t badRadius = 0, badAngle = 0, constrainedAngle = 0, tooLarge = 0;
            size_t freeBelowOne = 0, prescribedBelowOne = 0;
            double minAngle = 180.0, radiusEdge95 = 0.0, maxRadiusEdge = 0.0;

            bool targetsMet() const
            {
                return invalid == 0 && nearFlat == 0 && badRadius == 0 && badAngle == 0 && tooLarge == 0;
            }
        };

        struct PatchAttempts
        {
            std::array<int, 4> vertices = {};
            size_t visits = 0, accepted = 0, boundary = 0, region = 0, topology = 0, quality = 0, spacing = 0, placement = 0;
        };

        struct Work
        {
            size_t visits = 0, budgetStops = 0, offCenters = 0, orderedCandidates = 0, boundaryPoints = 0;
            const RTetraMeshProfile *profile = nullptr;
        };

        //! Export up to 50 worst cells, their neighbors, constraints and observed rejected attempts.
        static void writeReport(const RTetraMeshData &data, const RTetraMeshOptions &options,
                                const std::vector<PatchAttempts> &attempts, const std::string &path, const Work *work = nullptr);
        static Element measure(const double *a, const double *b, const double *c, const double *d, const bool *prescribedFaces = nullptr);
        //! Deterministic reduction; large meshes use up to four read-only workers.
        //! workers=1 forces serial evaluation; workers=0 chooses automatically.
        static Summary summarize(const RTetraMeshData &data, const RTetraMeshOptions &options, unsigned int workers = 0);
};
#endif
