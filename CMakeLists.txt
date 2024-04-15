cmake_minimum_required(VERSION 3.22)

project(lpn)

message (STATUS "Source dir ${PROJECT_SOURCE_DIR}")

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set (CMAKE_CXX_STANDARD 20)
    set (CMAKE_CXX_EXTENSIONS ON)
    set (CMAKE_CXX_STANDARD_REQUIRED ON)

    set (DEBUG_INFO_FLAGS "-g -gdwarf-4")
    set (CMAKE_CXX_FLAGS                     "${CMAKE_CXX_FLAGS} ${COMPILER_FLAGS}")
    set (CMAKE_CXX_FLAGS_RELWITHDEBINFO      "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O3 ${DEBUG_INFO_FLAGS} ${CMAKE_CXX_FLAGS_ADD}")
    set (CMAKE_CXX_FLAGS_DEBUG               "${CMAKE_CXX_FLAGS_DEBUG} -O0 ${DEBUG_INFO_FLAGS} ${CMAKE_CXX_FLAGS_ADD}")
endif()

set (CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(cmake/glob_sources.cmake)
include(cmake/sanitize.cmake)

include(CTest)

add_subdirectory(src)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Release")
    set (ASAN true)
    set (UBSAN true)
    add_subdirectory(tests)
endif()
add_subdirectory(programs)
