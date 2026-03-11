add_executable(test_mixed test_mixed_scalar.cpp test_mixed_binop.cpp test_mixed_unop.cpp
        test_mixed_comparison.cpp test_mixed_misc.cpp)
target_link_libraries(test_mixed caffeine winterval dualdomain GTest::gtest_main)
include(GoogleTest)
