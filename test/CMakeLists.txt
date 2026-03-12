include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# Windows compatibility
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../../../out/tests)

add_subdirectory(affine)
add_subdirectory(winterval)
add_subdirectory(mixedform)
add_subdirectory(dualnumber)
