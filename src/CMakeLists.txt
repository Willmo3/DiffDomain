set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../../out/)
add_subdirectory(Winterval)
add_subdirectory(Caffeine)
add_subdirectory(MixedForm)
add_subdirectory(SampleRange)

# Uncomment for separate real numbered library
# add_library(real Real/Real.cpp)

# Uncomment for unified domain library
add_library(domains Real/Real.cpp)
target_link_libraries(domains winterval caffeine mixedform samplerange)

