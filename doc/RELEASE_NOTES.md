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

- Quadratic (second order) elements can be generated: 10-node tetrahedra, 6-node
  triangles and 3-node edges with a shared node in the middle of each straight edge
- Volumes separated by interfaces and volumes nested in each other are meshed as
  separate regions
- Points and lines are kept, lines are split where nodes are added on them
- An existing mesh can be reconstructed: its volumes keep their markers and their
  interfaces are kept also if there are no triangles on them
- Mesh size can be prescribed at input nodes, it is interpolated to added nodes
- Triangles lying in one plane are treated as one facet, so the generator is not
  bound to the input triangulation of flat faces
- Mesh generation reports phase text through `RProgress`, including boundary
  recovery, refinement, optimization, output preparation and diagnostics. Refinement
  scans and optimization sweeps report percentage progress; adaptive queues pulse
  with work counts. Updates are throttled, status text is mirrored to the log without
  percentages. Phases with additional log output have a timed, indented section;
  phases without it print only their status line. Work-count updates stay within
  the current section. Progress and log indentation are finalized on success or
  exception
- **RBoundaryRecovery** can adopt a bounded planar patch's triangulation from the
  tetrahedral mesh before adding more vertices. The replacement is committed only
  after checking the patch boundary, constrained edges, vertex set, area, orientation
  and surface marker, allowing a valid planar patch to use the triangulation
  already present in the volume mesh
- Boundary recovery retries a failed circumcenter placement at the triangle centroid
  when placement is blocked by an existing vertex, an unusable facet location or an
  invalid cavity. Quality refinement enforces spacing and size safeguards
- Boundary recovery errors identify the actual split rejection, face or edge and node
  identifiers, marker, full-precision coordinates, circumradius or edge length, and
  minimum split length. Diagnostics distinguish recovery failures from detected
  surface intersections
