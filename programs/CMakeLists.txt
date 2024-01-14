find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_headers_and_sources(lpn ${CMAKE_CURRENT_SOURCE_DIR})
add_executable(lpn ${lpn_sources})
target_link_libraries(lpn PRIVATE lpn_core)
target_include_directories(lpn PRIVATE "${PROJECT_SOURCE_DIR}/src")
