find_package(Boost REQUIRED)

add_sources(lpn_thread ${CMAKE_CURRENT_SOURCE_DIR})
add_library(lpn_thread ${lpn_thread_sources})
