cmake_minimum_required(VERSION 3.5)
project(DiffDomain)
include_directories(lib)
include_directories(src)

set(CMAKE_CXX_STANDARD 20)

# Link-time optimization, so small accessors defined in .cpp files can inline across the domain libraries.
# CMP0069 must be NEW for CMake to honour IPO on non-Intel compilers given our old minimum version.
cmake_policy(SET CMP0069 NEW)
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported)
if (ipo_supported)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
endif()

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(lib)
