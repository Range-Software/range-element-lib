#ifndef REL_TETRA_MESH_REFINER_H
#define REL_TETRA_MESH_REFINER_H

#include <cstdint>

#include "rel_boundary_recovery.h"
#include "rel_tetra_mesh.h"
#include "rel_tetra_mesh_options.h"
#include "rel_tetra_mesh_profile.h"

//! Delaunay refinement of a tetrahedral mesh with recovered boundary.
//!
//! Tetrahedra which are badly shaped or too large are removed by inserting a vertex to the center
//! of their circumscribed sphere. If the center lies behind the boundary or too close to it,
//! the subsegment or the subface in the way is split instead.
class RTetraMeshRefiner
{

    protected:

        //! Quality of a tetrahedron.
        struct Quality
        {
            //! Tetrahedron has to be refined.
            bool bad = false;
            //! Tetrahedron is too large (otherwise it is badly shaped).
            bool tooLarge = false;
            //! Priority of refinement.
            double priority = 0.0;
            //! Center of circumscribed sphere.
            double center[3] = {0.0, 0.0, 0.0};
            //! Radius of circumscribed sphere.
            double radius = 0.0;
            //! Length of the shortest edge.
            double minEdge = 0.0;
            //! Length of edges the tetrahedron should have if it is too large.
            double targetEdge = 0.0;
        };

        //! Tetrahedron waiting to be refined.
        struct BadTetra
        {
            double priority;
            int tetraID;
            int v[4];
            int nAttempts;

            bool operator<(const BadTetra &badTetra) const
            {
                return this->priority < badTetra.priority;
            }
        };

        RTetraMesh &mesh;
        RBoundaryRecovery &boundary;
        const RTetraMeshOptions &options;

        //! Length below which nothing is split.
        double minLength;
        //! Number of added vertices at which refinement stops.
        uint64_t maxSteinerVertices;
        //! Whether the refinement was stopped before all tetrahedra were refined.
        bool stopped;
        uint64_t spacingRejections = 0;
        uint64_t boundaryRejections = 0;
        RTetraMeshProfile profile;

    private:

        //! Find quality of a tetrahedron.
        Quality findQuality(int tetraID) const;

        //! Refine tetrahedron. Return true if it should be tried again.
        bool refineTetra(const BadTetra &badTetra, const Quality &quality);

    public:

        //! Constructor.
        RTetraMeshRefiner(RTetraMesh &mesh, RBoundaryRecovery &boundary, const RTetraMeshOptions &options);

        //! Set length below which nothing is split.
        void setMinLength(double minLength);

        //! Set number of added vertices at which refinement stops.
        void setMaxSteinerVertices(uint64_t maxSteinerVertices);

        //! Return whether the refinement was stopped before all tetrahedra were refined.
        bool getStopped() const;

        uint64_t getSpacingRejections() const
        {
            return spacingRejections;
        }

        uint64_t getBoundaryRejections() const
        {
            return boundaryRejections;
        }

        const RTetraMeshProfile &getProfile() const
        {
            return profile;
        }

        //! Refine mesh.
        void refine();
};

#endif // REL_TETRA_MESH_REFINER_H
