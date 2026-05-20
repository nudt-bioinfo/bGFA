/*
 * @file /bGFA/src/bgfa_subcommand.hpp
 * @brief
 * @details
 * @author grampart
 * @version
 * @copyright: Copyright(c) 2025 by grampart, All Rights Reserved.
 * @license: MIT License
 */
#ifndef BGFA_SUBCOMMAND_HPP_
#define BGFA_SUBCOMMAND_HPP_

#include "bgfa_args.hpp"
#include "bgfa_graph.hpp"
#include "bgfa_seed.hpp"
#include "bgfa_merge.hpp"
#include <cctype>
#include <random>
#include <limits>
#include <fstream>

double ConvertProcess(CommandParser parse)
{
    struct timeval start_time, end_time;

    GFA gfa;
    // Determine input suffix and validate
    const std::string &in_path = parse.convert_args.gfa_path;
    size_t last_dot = in_path.find_last_of('.');
    if (last_dot == std::string::npos || last_dot == in_path.size() - 1)
    {
        std::cout << "[Convert] Input file has no valid extension: "
                  << in_path << std::endl;
    }

    std::string in_ext = in_path.substr(last_dot + 1);
    for (auto &ch : in_ext)
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));

    if (in_ext != "gfa" && in_ext != "bgfa")
    {
        std::cout << "[Convert] Unsupported input format: ."
                  << in_ext << ", expected .gfa or .bgfa" << std::endl;
    }

    // Auto output naming: gfa -> bgfa, bgfa -> gfa
    if (parse.convert_args.output_path == "auto")
    {
        const std::string base = in_path.substr(0, last_dot);
        parse.convert_args.output_path = (in_ext == "gfa") ? (base + ".bgfa") : (base + ".gfa");
    }
    std::cout << "Converting " << in_path << " to " << parse.convert_args.output_path << std::endl;

    if (parse.convert_args.path_mode.empty())
        parse.convert_args.path_mode = "inc";
    for (auto &ch : parse.convert_args.path_mode)
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));

    std::cout << "\t Path Mode: " << parse.convert_args.path_mode << std::endl;
    std::cout << "\t Segment ID Storage: "
              << (parse.convert_args.segment_no_id ? "disabled (order-assigned)" : "enabled")
              << std::endl;

    gettimeofday(&start_time, NULL);
    gfa = GFA(in_path, parse.convert_args.path_mode);
    if (in_ext == "gfa" && parse.convert_args.normalized)
    {
        uword_t merged = gfa.mergeUniqueSuccessorByReverseTable();
        std::cout << "\t Normalized merge applied, merged pairs: " << merged << std::endl;
    }

    // Determine conversion mode: default "all"
    std::string mode = parse.convert_args.output_mode;
    for (auto &ch : mode)
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    if (mode.empty())
        mode = "all";

    if (in_ext == "gfa")
    {
        if (mode == "split")
        {
            // Determine base path without extension from output or input
            std::string base_out;
            if (parse.convert_args.output_path == "auto")
            {
                base_out = in_path.substr(0, last_dot);
            }
            else
            {
                const std::string &out = parse.convert_args.output_path;
                size_t out_last_dot = out.find_last_of('.');
                if (out_last_dot != std::string::npos)
                {
                    std::string out_ext = out.substr(out_last_dot + 1);
                    for (auto &ch2 : out_ext)
                        ch2 = static_cast<char>(std::tolower(static_cast<unsigned char>(ch2)));
                    base_out = (out_ext == "bgfa") ? out.substr(0, out_last_dot) : out;
                }
                else
                {
                    base_out = out;
                }
            }

            // Write split files and main BGFA via GFA method
            gfa.write2BgfaSplit(base_out, !parse.convert_args.segment_no_id);
            // Ensure convert output_path reflects main bgfa
            parse.convert_args.output_path = base_out + ".bgfa";
        }
        else
        {
            // default: full conversion
            gfa.write2Bgfa(parse.convert_args.output_path, !parse.convert_args.segment_no_id);
        }
    }
    else
    {
        // bgfa -> gfa
        gfa.write2GFA(parse.convert_args.output_path);
    }
    gettimeofday(&end_time, NULL);

    double total_time = (double)(end_time.tv_sec - start_time.tv_sec) +
                        (double)(end_time.tv_usec - start_time.tv_usec) / 1e6;

    std::cout << "Data successfully written to " << parse.convert_args.output_path << std::endl;
    return total_time;
}

double IndexProcess(CommandParser parse)
{
    GFA gfa;
    struct timeval start_time, end_time;
    std::cout << "Indexing " << parse.index_args.input_path
              << ", [mode] " << parse.index_args.seed_mode
              << " [iteration mode]" << parse.index_args.iterate_mode
              << std::endl;

    if (parse.index_args.output_path == "auto")
    {
        size_t last_dot = parse.index_args.input_path.find_last_of('.');
        if (last_dot != std::string::npos)
            parse.index_args.output_path = parse.index_args.input_path.substr(0, last_dot) + ".bmin";
        else
            parse.index_args.output_path = parse.index_args.input_path + ".bmin";
    }
    std::cout << "Output path: " << parse.index_args.output_path << std::endl;

    gettimeofday(&start_time, NULL);

    gfa = GFA(parse.index_args.input_path);

    // TODO: 是否需要增加例如FM-index之类的索引结构
    if (parse.index_args.seed_mode == "kmer")
    {
        // TODO: Kmer种子未实现
    }
    else if (parse.index_args.seed_mode == "minimizer")
    {
        std::cout << "\tMinimizer length: " << parse.index_args.minimizer_length
                  << ", Window size: " << parse.index_args.window_size << std::endl;
        MinimizerSet ms(gfa,
                        parse.index_args.minimizer_length,
                        parse.index_args.window_size,
                        parse.index_args.iterate_mode);

        ms.writeSeed2bmin(parse.index_args.output_path);
    }
    gettimeofday(&end_time, NULL);
    std::cout << "Seed data successfully written to " << parse.index_args.output_path << std::endl;

    double total_time = (double)(end_time.tv_sec - start_time.tv_sec) +
                        (double)(end_time.tv_usec - start_time.tv_usec) / 1e6;
    return total_time;
}

double ViewProcess(CommandParser parse)
{
    struct timeval start_time, end_time;
    GFA gfa;
    gettimeofday(&start_time, NULL);
    if (parse.view_args.graph_path != "")
    {
        gfa = GFA(parse.view_args.graph_path);
        std::cout << gfa;
    }
    else if (parse.view_args.minimizer_path != "")
    {
        MinimizerSet ms;
        ms = MinimizerSet(parse.view_args.minimizer_path);
        std::cout << ms;
    }
    gettimeofday(&end_time, NULL);
    double total_time = (double)(end_time.tv_sec - start_time.tv_sec) +
                        (double)(end_time.tv_usec - start_time.tv_usec) / 1e6;
    return total_time;
}

double MergeProcess(CommandParser parse)
{
    struct timeval start_time, end_time;
    if (WORD_BIT < 32)
    {
        std::cerr << "Error: merge bGFA requires WORD_BIT >= 32." << std::endl;
        return 0;
    }

    std::vector<GFA> gfa_list;
    GFA merged_gfa;

    gettimeofday(&start_time, NULL);

    for (const auto &graph_path : parse.merge_args.graph_path)
    {
        GFA gfa(graph_path);
        if (parse.merge_args.merge_type == 1 &&
            parse.merge_args.details_type == 1)
            gfa.setSegmentsDis(); // set distance for distance-aided merge

        gfa_list.push_back(gfa);
    }

    if (gfa_list.empty())
        throw std::invalid_argument("No input graphs provided.");

    GFA ref = gfa_list[0];
    std::vector<GFA> others;
    for (size_t i = 1; i < gfa_list.size(); ++i)
        others.push_back(gfa_list[i]);

    if (parse.merge_args.merge_type == 0) // strict merge
    {
        std::cout << "Merging " << parse.merge_args.graph_path.size()
                  << " files into " << parse.merge_args.output_path
                  << " with strict match." << std::endl;
        // Use first graph as reference, merge others into it using exact merge

        mergeGFAStrict(ref, others);
    }
    else if (parse.merge_args.merge_type == 1)
    {
        std::cout << "Merging " << parse.merge_args.graph_path.size()
                  << " files into " << parse.merge_args.output_path;

        if (parse.merge_args.details_type == 0)
        {
            std::cout << " with exact match." << std::endl;

            mergeGFAExact(ref, others, parse.merge_args.merge_limited, parse.merge_args.alignment_kernel);
        }
        else if (parse.merge_args.details_type == 1)
        {
            std::cout << " with distance-aided exact match." << std::endl;

            mergeGFAExactByDistance(ref, others, parse.merge_args.merge_limited, parse.merge_args.alignment_kernel);
        }
        else if (parse.merge_args.details_type == 2)
        {
            std::cout << " with minimizer-aided exact match." << std::endl;

            mergeGFAExactByMinimize(ref, others, parse.merge_args.merge_limited, parse.merge_args.alignment_kernel);
        }
        else
        {
            throw std::invalid_argument("Invalid details type: " + std::to_string(parse.merge_args.details_type));
        }
    }
    else
    {
        throw std::invalid_argument("Invalid merge type: " + std::to_string(parse.merge_args.merge_type));
    }

    ref.write2Bgfa(parse.merge_args.output_path);

    gettimeofday(&end_time, NULL);

    double total_time = (double)(end_time.tv_sec - start_time.tv_sec) +
                        (double)(end_time.tv_usec - start_time.tv_usec) / 1e6;
    return total_time;
}

double StatsProcess(CommandParser parse)
{
    struct timeval start_time, end_time;
    gettimeofday(&start_time, NULL);

    GFA gfa(parse.view_args.graph_path);
    gfa.briefStat();

    gettimeofday(&end_time, NULL);
    double total_time = (double)(end_time.tv_sec - start_time.tv_sec) +
                        (double)(end_time.tv_usec - start_time.tv_usec) / 1e6;
    return total_time;
}

#endif // BGFA_SUBCOMMAND_HPP_