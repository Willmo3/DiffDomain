# NOTE: due to static state, ordering of these test cases is important!
add_executable(test_caffeine test_affine_binop.cpp
        test_affine_unop.cpp
        test_affine_scalar.cpp
        test_affine_comparison.cpp
        test_affine_misc.cpp)
target_link_libraries(test_caffeine caffeine winterval GTest::gtest_main)
include(GoogleTest)

