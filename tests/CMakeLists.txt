# Unit tests for range-element-lib (QTest based).
#
# Enabled by configuring the project with -DRANGE_BUILD_TESTS=ON.
# Run with `ctest` from the build directory.

find_package(Qt6 REQUIRED COMPONENTS Test)

set(CMAKE_AUTOMOC ON)

set(range_element_lib_tests
    tst_geometric_predicates
    tst_delaunay_tetrahedralizer
    tst_tetra_mesher
    tst_tetra_mesh_optimizer
)

foreach(test_name IN LISTS range_element_lib_tests)
    qt_add_executable(${test_name} ${test_name}.cpp)
    target_link_libraries(${test_name}
        PRIVATE
            range-element-lib
            Qt6::Test
    )
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()

# Fail non-converging refinement regressions instead of leaving CI running.
set_tests_properties(tst_tetra_mesher PROPERTIES TIMEOUT 60)
