# Range Element Library

`range-element-lib` generates finite elements. Its main purpose is to generate a 3D tetrahedral mesh of volumes enclosed by a triangulated surface, which `range-model-lib` uses to mesh an `RModel`.

The library does not know `RModel`. It works with a plain description of a mesh (`RTetraMeshData`), so that it depends only on `range-base-lib`. Conversion between `RModel` and `RTetraMeshData` is done by `RModelTetraMesh` in `range-model-lib`.

## How a mesh is generated

`RTetraMesher::generate()` runs these steps:

1. **Delaunay tetrahedralization** of all input nodes (`RDelaunayTetrahedralizer`, Bowyer-Watson algorithm).
2. **Boundary recovery** (`RBoundaryRecovery`). Input edges and triangles which are missing in the tetrahedralization are recovered by adding nodes on them.
3. **Removal of tetrahedra outside of the domain**. Each group of tetrahedra enclosed by triangles becomes a volume (region), tetrahedra which can be reached from outside are removed.
4. **Refinement** (`RTetraMeshRefiner`, optional). Tetrahedra which are badly shaped or too large are removed by adding a node to the center of their circumscribed sphere, or by splitting the boundary which is in the way.

## Classes

| Class | Role |
|---|---|
| `RTetraMesher` | Entry point, `generate(input, options, output)` |
| `RTetraMeshData` | Input and output: nodes, edges, triangles, tetrahedra and their markers |
| `RTetraMeshOptions` | Options: quality mesh, radius/edge ratio, minimum dihedral angle, maximum volume, size function, tolerance, limit of added nodes |
| `RGeometricPredicates` | Exact `orient3d()` and `inSphere()` tests |
| `RTetraMesh` | Tetrahedralization which points can be inserted to: point location, cavity search, cavity filling |
| `RDelaunayTetrahedralizer` | Delaunay tetrahedralization of a set of points |
| `RBoundaryRecovery` | Subsegments and subfaces of the boundary, their recovery and splitting |
| `RTetraMeshRefiner` | Refinement by quality and size criteria |

## Input

- **Nodes** with coordinates, an optional marker of a point and an optional required mesh size.
- **Triangles** of the boundary and of interfaces between volumes, with a marker of a surface.
- **Edges** which must be present in the mesh, with a marker of a line.
- **Tetrahedra** of an existing mesh with markers of volumes, used only if the mesh is reconstructed. Their faces between different volumes become interfaces, also if there are no triangles on them.

The surface has to be **conforming**: triangles may touch only in whole edges or in nodes. Triangles which intersect each other, or a node of one triangle lying inside of an edge or inside of another triangle, are not supported and make the generation fail with an error. Nodes with the same coordinates are merged.

## Output

- Input nodes at the same positions, followed by added nodes.
- Triangles of all surfaces and interfaces. Triangles lying in one plane and having the same marker form a facet, which triangulation is free to change.
- Edges with a marker, split where nodes were added on them.
- Tetrahedra with a marker of their volume and with their neighbors. Node order gives a positive volume.

## Robustness

- Geometric tests are exact. They are evaluated in floating point arithmetic with an error bound and repeated with floating point expansions when the bound is not met (J. R. Shewchuk's method). Points lying on one sphere are handled by symbolic perturbation (O. Devillers, M. Teillaud).
- A new point is always inserted to a cavity which it can see the whole boundary of, so no tetrahedron with zero or negative volume is ever created.
- Once the boundary is recovered, a cavity never crosses it. Nodes added on one surface therefore can not damage another surface, however close it is.
- Points closer to a plane than the tolerance are treated as lying in it. Faces which are planar only up to rounding errors (a plane which is not aligned with axes) therefore do not produce flat tetrahedra.

## Limitations

- Refinement does not shorten edges where parts of the boundary meet in a sharp angle. Elements there may not meet the required radius/edge ratio or dihedral angle.
- Flat tetrahedra (slivers) are removed only by refinement (`RTetraMeshOptions::setMinDihedralAngle()`). Slivers touching the boundary may remain.
- Surfaces with errors (intersections, gaps, hanging nodes) are not repaired.

## Tests

Unit tests are built with `-DRANGE_BUILD_TESTS=ON` and run by `ctest`.

## Class prefix

All public classes use the `R` prefix, source files use the `rel_` prefix.

## Dependencies

- `range-base-lib`
- Qt 6.10+
- C++20
