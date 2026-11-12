## Version 1.0.0

### New library

`range-element-lib` generates a 3D tetrahedral mesh of volumes enclosed by a
triangulated surface.

- **RTetraMesher** generates the mesh: Delaunay tetrahedralization of input nodes,
  recovery of input edges and triangles, removal of tetrahedra outside of the
  domain and refinement
- **RTetraMeshData** holds the input (nodes, edges, triangles and tetrahedra of a
  mesh to be reconstructed) and the output of the generator, with markers of
  points, lines, surfaces and volumes
- **RTetraMeshOptions** sets quality mesh, maximum radius/edge ratio, minimum
  dihedral angle, maximum element volume, use of mesh sizes prescribed at nodes,
  tolerance of coplanarity and a limit of added nodes
- **RGeometricPredicates** provides exact orientation and in-sphere tests
- **RTetraMesh**, **RDelaunayTetrahedralizer**, **RBoundaryRecovery** and
  **RTetraMeshRefiner** are the steps of mesh generation, which can be used on
  their own

### Features

- Volumes separated by interfaces and volumes nested in each other are meshed as
  separate regions
- Points and lines are kept, lines are split where nodes are added on them
- An existing mesh can be reconstructed: its volumes keep their markers and their
  interfaces are kept also if there are no triangles on them
- Mesh size can be prescribed at input nodes, it is interpolated to added nodes
- Triangles lying in one plane are treated as one facet, so the generator is not
  bound to the input triangulation of flat faces
