#ifndef REL_TETRA_MESHER_H
#define REL_TETRA_MESHER_H

#include "rel_tetra_mesh_data.h"
#include "rel_tetra_mesh_options.h"
#include "rel_tetra_mesh_profile.h"

//! Tetrahedral mesh generator.
//!
//! Mesh is generated in following steps:
//! 1. Delaunay tetrahedralization of all input nodes.
//! 2. Recovery of input edges and triangles, which may add nodes on them.
//! 3. Removal of tetrahedra lying outside of volumes enclosed by triangles.
//! 4. Refinement of tetrahedra which are badly shaped or too large (optional).
class RTetraMesher
{

    private:

        //! Constructor.
        RTetraMesher() {}

    public:

        //! Generate tetrahedral mesh.
        //! Throws RError if input is not valid or if mesh could not be generated.
        static void generate(const RTetraMeshData &input, const RTetraMeshOptions &options, RTetraMeshData &output);
        //! Generate and optionally return phase timings and optimizer work counters.
        static void generate(const RTetraMeshData &input, const RTetraMeshOptions &options, RTetraMeshData &output, RTetraMeshProfile *profile);

};

#endif // REL_TETRA_MESHER_H
