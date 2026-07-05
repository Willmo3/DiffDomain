set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../../out/)
add_subdirectory(Winterval)
add_subdirectory(Caffeine)
add_subdirectory(MixedForm)
add_subdirectory(SampleRange)

# Uncomment for separate real numbered library
# add_library(real Real/Real.cpp)

# Uncomment for unified domain library
add_library(domains Real/Real.cpp)
# Note: we statically link polylib to simplify build process, but this may result in multiple libs bringing the package in
# Additionally, note that this requires polylib to be recognized before it can be built.
target_link_libraries(domains winterval caffeine mixedform samplerange polylib)

