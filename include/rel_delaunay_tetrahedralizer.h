#ifndef REL_DELAUNAY_TETRAHEDRALIZER_H
#define REL_DELAUNAY_TETRAHEDRALIZER_H

#include <vector>

#include "rel_tetra_mesh.h"

//! Builds Delaunay tetrahedralization of points by inserting them one by one (Bowyer-Watson algorithm).
class RDelaunayTetrahedralizer
{

    private:

        //! Constructor.
        RDelaunayTetrahedralizer()
        {
        }

    public:

        //! Tetrahedralize all points of a mesh which has no tetrahedra yet.
        //! Points with the same coordinates as an already inserted point, or too close to it to form
        //! tetrahedra with it, are not inserted. Position of the inserted point is set to duplicates
        //! for them, -1 is set for all other points.
        static void tetrahedralize(RTetraMesh &mesh, std::vector<int> &duplicates);
};

#endif // REL_DELAUNAY_TETRAHEDRALIZER_H
