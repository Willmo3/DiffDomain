add_executable(test_winterval test_winterval.cpp)
target_link_libraries(test_winterval winterval GTest::gtest_main)
include(GoogleTest)