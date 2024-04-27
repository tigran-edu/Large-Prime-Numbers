find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_sources(lpn_core ${CMAKE_CURRENT_SOURCE_DIR})
add_library(lpn_core ${lpn_core_sources})

target_include_directories(lpn_core PRIVATE "${PROJECT_SOURCE_DIR}/thread")
target_link_libraries(lpn_core PRIVATE lpn_thread)
