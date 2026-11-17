#ifndef REL_TETRA_MESH_STATISTICS_H
#define REL_TETRA_MESH_STATISTICS_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>

#include "rel_tetra_mesh_data.h"

//! Counts, quality ranges and TetGen-compatible histogram bins for a generated mesh.
//! Geometry is measured from corner nodes, also for quadratic elements.
class RTetraMeshStatistics
{
    public:

        struct Range
        {
            double minimum = std::numeric_limits<double>::infinity();
            double maximum = -std::numeric_limits<double>::infinity();

            void add(double value);
        };

        struct Statistics
        {
            size_t faces = 0;
            size_t exteriorFaces = 0;
            size_t validElements = 0;
            size_t invalidElements = 0;
            Range volume;
            Range edge;
            Range aspectRatio;
            Range faceAngle;
            Range dihedralAngle;
            std::array<uint64_t, 12> aspectHistogram = {};
            std::array<uint64_t, 18> faceAngleHistogram = {};
            std::array<uint64_t, 18> dihedralAngleHistogram = {};
            double minimumInputFaceAngle = std::numeric_limits<double>::infinity();
            double minimumInputDihedralAngle = std::numeric_limits<double>::infinity();
        };

        //! Counts retained by the generator which are not encoded in its output data.
        struct Counts
        {
            size_t inputSegments = 0;
            size_t inputRegions = 0;
            size_t meshSegments = 0;
            size_t meshRegions = 0;
            uint64_t steinerOnFacets = 0;
            uint64_t steinerOnSegments = 0;
            uint64_t steinerInside = 0;
            size_t midEdgeNodes = 0;
        };

        //! Measure valid, positively oriented tetrahedra. Invalid cells are counted
        //! separately and excluded from quality ranges and histograms.
        static Statistics measure(const RTetraMeshData &input, const RTetraMeshData &output);

        //! Print Range-style indented sections and aligned plain-text tables.
        //! Each row is one complete logger message.
        static void print(const RTetraMeshData &input, const RTetraMeshData &output,
                          const Counts &counts, const Statistics &statistics);
};

#endif
