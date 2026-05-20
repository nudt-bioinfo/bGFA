/*
 * @file /bGFA/src/bgfa_alignment.hpp
 * @brief
 * @details
 * @author grampart
 * @version
 * @copyright: Copyright(c) 2025 by grampart, All Rights Reserved.
 * @license: MIT License
 */
#ifndef BGFA_ALIGNMENT_HPP_
#define BGFA_ALIGNMENT_HPP_

#include "bgfa_graph.hpp"

#include "ksw/ksw2.h"
#include "wfa/bindings/cpp/WFAligner.hpp"
using namespace wfa;

/**
 * @brief Convert KSW CIGAR string to MX format
 * @param {string} &cigar_str_in - Input CIGAR string from KSW
 * @param {string} &tseq - Target sequence
 * @param {string} &qseq - Query sequence
 * @return {string} Converted CIGAR string in MX format
 */
std::string convertKswCigar2MX(const std::string &cigar_str_in,
                               const std::string &tseq, const std::string &qseq)
{
    if (cigar_str_in.empty())
        return cigar_str_in;

    std::string new_cigar;
    new_cigar.reserve(cigar_str_in.size());

    size_t p = 0, n = cigar_str_in.size();
    size_t pos_t = 0, pos_q = 0; // positions in tseq / qseq

    auto append_op = [&](uint64_t len, char op)
    {
        if (len == 0)
            return;
        if (!new_cigar.empty())
        {
            // if previous op char equals current op, merge lengths
            char last_op = new_cigar.back();
            if (last_op == op)
            {
                // find start index of last token's number (the digits before last_op)
                size_t pos_op = new_cigar.size() - 1;
                size_t num_start = pos_op;
                while (num_start > 0 && std::isdigit(static_cast<unsigned char>(new_cigar[num_start - 1])))
                    --num_start;
                size_t num_len = pos_op - num_start;
                uint64_t last_len = 0;
                if (num_len > 0)
                {
                    last_len = std::stoull(new_cigar.substr(num_start, num_len));
                }
                uint64_t merged = last_len + len;
                new_cigar.erase(num_start);
                new_cigar += std::to_string(merged);
                new_cigar.push_back(op);
                return;
            }
        }
        new_cigar += std::to_string(len);
        new_cigar.push_back(op);
    };

    while (p < n)
    {
        uint64_t num = 0;
        while (p < n && std::isdigit(static_cast<unsigned char>(cigar_str_in[p])))
        {
            num = num * 10 + (cigar_str_in[p] - '0');
            ++p;
        }
        if (p >= n)
            break; // malformed
        char op = cigar_str_in[p++];

        if (op == 'M')
        {
            // split M into runs of matches (M) and mismatches (X)
            uint64_t run_len = 0;
            bool run_match = false;
            for (uint64_t i = 0; i < num; ++i)
            {
                char a = (pos_t + i < tseq.size()) ? tseq[pos_t + i] : '\0';
                char b = (pos_q + i < qseq.size()) ? qseq[pos_q + i] : '\0';
                bool equal = (std::toupper(static_cast<unsigned char>(a)) ==
                              std::toupper(static_cast<unsigned char>(b)));
                if (i == 0)
                {
                    run_match = equal;
                    run_len = 1;
                }
                else if (equal == run_match)
                {
                    ++run_len;
                }
                else
                {
                    append_op(run_len, run_match ? 'M' : 'X');
                    run_match = equal;
                    run_len = 1;
                }
            }
            if (run_len)
                append_op(run_len, run_match ? 'M' : 'X');

            pos_t += num;
            pos_q += num;
        }
        else if (op == 'I')
        {
            append_op(num, 'I');
            pos_q += num;
        }
        else if (op == 'D')
        {
            append_op(num, 'D');
            pos_t += num;
        }
        else
        {
            // keep other ops as-is
            append_op(num, op);
        }
    }

    return new_cigar;
}

/**
 * @brief: Use KSW2 to perform sequence alignment and return the CIGAR string
 * @param {string} &tseq - Target sequence
 * @param {string} &qseq - Query sequence
 * @return {string} CIGAR string
 */
std::string alignmentSequenceKSW(const std::string &tseq,
                                 const std::string &qseq)
{

    std::string cigar_str;
    if (tseq.empty() || qseq.empty())
    {
        return "0M";
    }

    int i, a = 2, b = 4; // a>0 and b<0
    int8_t mat[25] = {a, b, b, b, 0,
                      b, a, b, b, 0,
                      b, b, a, b, 0,
                      b, b, b, a, 0,
                      0, 0, 0, 0, 0};
    int tl = tseq.size(), ql = qseq.size();
    uint8_t *ts, *qs, c[256];
    ksw_extz_t ez;

    memset(&ez, 0, sizeof(ksw_extz_t));
    memset(c, 4, 256);
    c['A'] = c['a'] = 0;
    c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2;
    c['T'] = c['t'] = 3; // build the encoding table
    ts = (uint8_t *)malloc(tl);
    qs = (uint8_t *)malloc(ql);
    for (i = 0; i < tl; ++i)
        ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
    for (i = 0; i < ql; ++i)
        qs[i] = c[(uint8_t)qseq[i]];

    ksw_extz2_sse(0, ql, qs, tl, ts, 5, mat, 4, 2, -1, -1, 0, 0, &ez);

    if (ez.n_cigar == 0)
    {
        free(ts);
        free(qs);
        return "0M";
    }
    cigar_str.reserve(ez.n_cigar * 4); // Reserve space for CIGAR string
    for (i = 0; i < ez.n_cigar; ++i)   // convert CIGAR to string
    {
        int len = ez.cigar[i] >> 4;
        char op = GlobalVariant::bin2cigarop[ez.cigar[i] & 0xf];
        cigar_str += std::to_string(len) + op;
    }
    free(ez.cigar);
    free(ts);
    free(qs);

    // Convert KSW CIGAR (M include mismatches) into M/X-aware CIGAR
    cigar_str = convertKswCigar2MX(cigar_str, tseq, qseq);

#ifdef INFO_DEBUG
    std::cout << "[KSW Alignment] CIGAR: " << cigar_str
              << "\n\t\tTarget Seq: " << tseq
              << "\n\t\tQuery  Seq: " << qseq << std::endl;
#endif

    return cigar_str;
}

/**
 * @brief: Use WFA to perform sequence alignment and return the CIGAR string
 * @param {string} &tseq - Target sequence
 * @param {string} &qseq - Query sequence
 * @return {string} CIGAR string
 */
std::string alignmentSequenceWFA(const std::string &tseq,
                                 const std::string &qseq)
{

    std::string cigar_str;
    if (tseq.empty() || qseq.empty())
    {
        return "0M"; // Return empty CIGAR string if either sequence is empty
    }

    WFAlignerGapAffine aligner(4, 6, 2, WFAligner::Alignment,
                               WFAligner::MemoryHigh);

    // Normalize sequences to uppercase to ensure consistent behavior
    // across aligner backends (WFA expects canonical bases).
    std::string t_up = tseq;
    std::string q_up = qseq;
    for (char &c : t_up)
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    for (char &c : q_up)
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));

    aligner.alignEnd2End(t_up, q_up);

    cigar_str = aligner.getCIGAR(true);
    std::replace(cigar_str.begin(), cigar_str.end(), '=', 'M');

#ifdef INFO_DEBUG
    std::cout << "[WFA Alignment] CIGAR:      " << cigar_str
              << "\n\t\tTarget Seq: " << tseq
              << "\n\t\tQuery Seq:  " << qseq << std::endl;
#endif

    return cigar_str;
}

std::string alignmentSequence(const std::string &tseq,
                              const std::string &qseq,
                              const uint8_t kernel = 1)
{
    if (kernel == 0) // Use KSW alignment kernel
    {
        return alignmentSequenceKSW(tseq, qseq);
    }
    else if (kernel == 1) // Use WFA alignment kernel
    {
        return alignmentSequenceWFA(tseq, qseq);
    }
    else
    {
        return "0M";
    }
}

#endif // BGFA_ALIGNMNET_HPP_