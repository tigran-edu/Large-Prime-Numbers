find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_executable(lpn main.cpp)
target_link_libraries(lpn PRIVATE lpn_core)
target_include_directories(lpn PRIVATE "${PROJECT_SOURCE_DIR}/src")
