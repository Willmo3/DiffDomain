add_executable(test_samplerange
        test_samplerange.cpp)
target_link_libraries(test_samplerange samplerange GTest::gtest_main)
include(GoogleTest)

