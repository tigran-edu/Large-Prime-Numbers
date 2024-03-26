find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_sources(lpn_core ${CMAKE_CURRENT_SOURCE_DIR})
add_library(lpn_core ${lpn_core_sources})
