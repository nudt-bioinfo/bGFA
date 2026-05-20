/*
 * @file /bGFA/src/test/test_bgfa_merge.hpp
 * @brief
 * @details
 * @author grampart
 * @version
 * @copyright: Copyright(c) 2025 by grampart, All Rights Reserved.
 * @license: MIT License
 */
#pragma once
#include "test_bgfa.hpp"

static std::string mergeDir = "./data/graph/merge";
static std::string s1 = "CTTCTTCAGATTACCCGCGTTGTAGTGCGATTATCTTGAGGTCTATATAGGTCGATCTATTATTTATGGTACTGAACCTCACCCAGACCCCCCAACCAAT";
static std::string s2 = "CTTCTTCAGAGTACCTGCGTTGTAGTGCGATTCTCTTGGGTCTATATAGGTCGATCTATTATTTTATGGTACTGAACCTCACCCAGCACCCCCCCAACCAAT";

static int test_alignment_kernel_ksw()
{
    std::string cigar = alignmentSequence(s1, s2, /*kernel=*/0); // KSW
    if (cigar.empty() || cigar == "0M")
    {
        std::cout << "[Failed] alignmentSequence returned empty/0M\n";
        return 0;
    }
    else
    {
#ifdef INFO_DEBUG
        std::cout << "KSW CIGAR: " << cigar << std::endl;
#endif
        std::cout << "\t[PASS]";
        return 1;
    }
}

static int test_alignment_kernel_wfa()
{
    std::string cigar = alignmentSequence(s1, s2, /*kernel=*/1); // WFA
    if (cigar.empty() || cigar == "0M")
    {
        std::cout << "[Failed] alignmentSequence returned empty/0M\n";
        return 0;
    }
    else
    {
#ifdef INFO_DEBUG
        std::cout << "WFA CIGAR: " << cigar << std::endl;
#endif
        std::cout << "\t[PASS]";
        return 1;
    }
}

static int test_merge_strict()
{
    GFA res(mergeDir + "/r3.gfa");
    GFA ref(mergeDir + "/g1.gfa");

    std::vector<GFA> others;
    others.emplace_back(GFA(mergeDir + "/g2.gfa"));

    mergeGFAStrict(ref, others);

    ref.compactAndRenumber();

#ifdef INFO_DEBUG
    std::string outGfa = mergeDir + "/res3.gfa";
    ref.write2GFA(outGfa);
#endif

    if (res == ref)
    {
        std::cout << "\t[PASS]";
        return 1;
    }
    else
    {
        std::cout << "[Failed] exact-merge structural check failed\n";
        return 0;
    }
}

static int test_merge_exact()
{
    GFA res(mergeDir + "/r4.gfa");
    GFA ref(mergeDir + "/g1.gfa");
    std::vector<GFA> others;
    others.emplace_back(GFA(mergeDir + "/g2.gfa"));

    mergeGFAExact(ref, others, /*merge_limited=*/90, /*alignment_kernel=*/1);

    ref.compactAndRenumber();

#ifdef INFO_DEBUG
    std::string outGfa = mergeDir + "/res4.gfa";
    ref.write2GFA(outGfa);
#endif

    if (res == ref)
    {
        std::cout << "\t[PASS]";
        return 1;
    }
    else
    {
        std::cout << "[Failed] exact-merge structural check failed\n";
        return 0;
    }
}

static int test_merge_strict_dis()
{
    GFA res(mergeDir + "/r5.gfa");
    GFA ref(mergeDir + "/g3.gfa");
    std::vector<GFA> others;
    others.emplace_back(GFA(mergeDir + "/g4.gfa"));

    mergeGFAStrictByDistance(ref, others);

    ref.compactAndRenumber();

#ifdef INFO_DEBUG
    std::string outGfa = mergeDir + "/res5.gfa";
    ref.write2GFA(outGfa);
#endif

    if (res == ref)
    {
        std::cout << "\t[PASS]";
        return 1;
    }
    else
    {
        std::cout
            << "[Failed] strict-merge-by-distance structural check failed\n";
        return 0;
    }
}

static int test_merge_exact_dis()
{
    // GFA res(mergeDir + "/r6.gfa");
    GFA ref(mergeDir + "/g3.gfa");
    std::vector<GFA> others;
    others.emplace_back(GFA(mergeDir + "/g4.gfa"));

    mergeGFAExactByDistance(ref, others, 90, 1);

    std::cout << ref;

    ref.compactAndRenumber();

#ifdef INFO_DEBUG
    std::string outGfa = mergeDir + "/res6.gfa";
    ref.write2GFA(outGfa);
#endif

    // if (res == ref)
    // {
    //     std::cout << "\t[PASS]";
    //     return 1;
    // }
    // else
    // {
    std::cout
        << "[Failed] strict-merge-by-distance structural check failed\n";
    return 0;
    // }
}

inline void test_bgfa_merge(int &pass_num, int &test_num)
{
    test_num += 6;

    std::cout << "Running BGFA Merge Tests..." << std::endl;

    pass_num += test_alignment_kernel_ksw();
    std::cout << " KSW Tests" << std::endl;

    pass_num += test_alignment_kernel_wfa();
    std::cout << " WFA Tests" << std::endl;

    pass_num += test_merge_strict();
    std::cout << " STRICT Merge" << std::endl;

    pass_num += test_merge_exact();
    std::cout << " EXACT Merge" << std::endl;

    pass_num += test_merge_strict_dis();
    std::cout << " STRICT Merge by Distance" << std::endl;

    pass_num += test_merge_exact_dis();
    std::cout << " EXACT Merge by Distance" << std::endl;

    // pass_num
}
