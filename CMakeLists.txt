# Create a static library
qt_add_library(range-element-lib
    STATIC
        src/rel_boundary_recovery.cpp
        src/rel_delaunay_tetrahedralizer.cpp
        src/rel_geometric_predicates.cpp
        src/rel_tetra_mesh.cpp
        src/rel_tetra_mesh_data.cpp
        src/rel_tetra_mesh_options.cpp
        src/rel_tetra_mesh_refiner.cpp
        src/rel_tetra_mesh_progress.h
        src/rel_tetra_mesher.cpp
        src/rel_tetra_mesh_statistics.cpp
        src/rel_tetra_mesh_quality.cpp
        src/rel_tetra_mesh_optimizer.cpp

        include/rel_boundary_recovery.h
        include/rel_delaunay_tetrahedralizer.h
        include/rel_geometric_predicates.h
        include/rel_tetra_mesh.h
        include/rel_tetra_mesh_data.h
        include/rel_tetra_mesh_options.h
        include/rel_tetra_mesh_refiner.h
        include/rel_tetra_mesher.h
        include/rel_tetra_mesh_statistics.h
        include/rel_tetra_mesh_quality.h
        include/rel_tetra_mesh_profile.h
        include/rel_tetra_mesh_optimizer.h
)

# Exact geometric predicates rely on IEEE 754 arithmetic: a product and a sum must not be
# contracted to a fused multiply-add and floating point expressions must not be reordered.
if(MSVC)
    target_compile_options(range-element-lib PRIVATE /fp:strict)
else()
    target_compile_options(range-element-lib PRIVATE -ffp-contract=off -fno-fast-math)
endif()

# Include directories
target_include_directories(range-element-lib
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(range-element-lib
    PUBLIC
        range-base-lib
        common_defines
)

# Unit tests (opt-in via -DRANGE_BUILD_TESTS=ON)
if(RANGE_BUILD_TESTS)
    add_subdirectory(tests)
endif()
