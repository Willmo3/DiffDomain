cmake_minimum_required(VERSION 3.5)
project(DiffDomain)
include_directories(lib)
include_directories(src)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(src)
add_subdirectory(test)