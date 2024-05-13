cmake_minimum_required(VERSION 3.22)

add_sources(lpn_unit_tests ${CMAKE_CURRENT_SOURCE_DIR})

find_package(GTest REQUIRED)

macro(add_own_test test_name)
    add_executable(${test_name} ${test_name}.cpp)
    add_test(${test_name} ${test_name})
    target_link_libraries(${test_name} PRIVATE lpn_core GTest::gtest GTest::gtest_main)
    target_include_directories(${test_name} PRIVATE "${PROJECT_SOURCE_DIR}/src")
endmacro()

add_own_test(basic_test)
add_own_test(rho_factor_test)
add_own_test(qs_factor_test)
add_own_test(gaussian_test)
add_own_test(congruence_test)
add_own_test(cfrac_test)
