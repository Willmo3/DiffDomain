add_executable(test_dualnumber
        test_dual_binop.cpp
        test_dual_scalar.cpp
        test_dual_unop.cpp
        test_dual_comparison.cpp
        test_dual_misc.cpp)
target_link_libraries(test_dualnumber winterval GTest::gtest_main)
include(GoogleTest)

