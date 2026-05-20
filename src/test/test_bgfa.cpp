/*
 * @file \bGFA\src\test\test_bgfa.cpp
 * @brief
 * @details
 * @author grampart
 * @version
 * @copyright: Copyright(c) 2025 by grampart, All Rights Reserved.
 * @license: MIT License
 */
#include "test_bgfa_basic.hpp"
#include "test_bgfa_merge.hpp"

int main()
{
    int pass_num = 0, test_num = 0;

    test_bgfa_basic(pass_num, test_num);
    // test_bgfa_merge(pass_num, test_num);

    std::cout << "BGFA Basic Graph Tests: Passed " << pass_num << " / "
              << test_num << std::endl;
}