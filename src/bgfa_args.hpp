/*
 * @Author: grampart guo018248@163.com
 * @Date: 2025-04-07 16:29:33
 * @LastEditors Please set LastEditors
 * @LastEditTime 2026-04-13 15:07:53
 * @FilePath: /bGFA/src/args.hpp
 */
#ifndef BGFA_ARGS_HPP_
#define BGFA_ARGS_HPP_

#include "./CLI11/CLI11.hpp"
#include <string>
#include <map>
#include <stdexcept>

typedef std::map<std::string, int> argsMap;
argsMap index_itrmode_map{{"single", 0}, {"single-reverse", 1}, {"double", 2}};
argsMap merge_alig_map{{"ksw", 0}, {"wfa", 1}};

struct ConvertArgs
{
    std::string gfa_path;
    std::string output_path = "auto";
    std::string output_mode = "all"; // "all" or "split"
    std::string path_mode = "inc";   // "direct" or "inc"
    bool normalized = false;
    bool segment_no_id = false;
    bool no_sosr = false;
    bool sosr_compact = false;
    // bool auto_output = true;
};

struct IndexArgs
{
    std::string input_path;
    std::string output_path;
    bool is_bgfa = false;
    std::string seed_mode;
    uint8_t iterate_mode;
    uint8_t window_size = 20;
    uint8_t minimizer_length = 15;
};

struct ViewArgs
{
    std::string graph_path;
    std::string minimizer_path;
};

struct MergeArgs
{
    std::vector<std::string> graph_path;
    uint8_t merge_type;           ///> 0: strict, 1: exactly
    uint8_t merge_limited = 90;   ///> percentage threshold for exact merge
    uint8_t alignment_kernel = 1; ///> 0: ksw, 1: wfa
    uint8_t details_type = 0;     ///> 0: none, 1: distance-aided, 2: minimizer-aided
    std::string output_path;
    // std::string minimizer_path;
};

class CommandParser
{
public:
    bool time = false;
    ConvertArgs convert_args;
    IndexArgs index_args;
    ViewArgs view_args;
    MergeArgs merge_args;
};

#endif // _ARGS_HPP_