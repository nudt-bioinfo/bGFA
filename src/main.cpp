/*
 * @Author: grampart guo018248@163.com
 * @Date: 2025-03-25 15:59:20
 * @LastEditors Please set LastEditors
 * @LastEditTime 2026-03-10 15:29:49
 * @FilePath: \bgfa\src\main.cpp
 */
#include "bgfa_subcommand.hpp"

int main(int argc, char *argv[])
{

#ifdef INFO_DEBUG
    std::cout << "DEBUG INFO is ON" << std::endl;
    std::cout << "THE BASIC SIZE IS " << WORD_BIT << std::endl;
#endif

    CommandParser parse;

    double time = 0;

    CLI::App app{"bGFA — a small genomics graph toolkit"};
    app.set_version_flag("-v, --version", "1.0.0");
    app.add_flag("-t, --time", parse.time, "Enable time info output");

    // convert
    auto &convert = *app.add_subcommand("convert", "Convert GFA to BGFA");
    convert.add_option("-g,--gfa", parse.convert_args.gfa_path,
                       "Input GFA file")
        ->required()
        ->check(CLI::ExistingFile);
    convert.add_option("-o,--output", parse.convert_args.output_path,
                       "Output BGFA file")
        ->default_str("auto");
    convert.add_option("--output-mode", parse.convert_args.output_mode,
                       "Conversion mode (all, split)")
        ->default_str("all")
        ->check(CLI::IsMember({"all", "split"}, CLI::ignore_case));
    convert.add_option("--path-mode", parse.convert_args.path_mode,
                       "Path encoding mode (direct, inc)")
        ->default_str("inc")
        ->check(CLI::IsMember({"direct", "inc"}, CLI::ignore_case));
    convert.add_flag("--normalized", parse.convert_args.normalized,
                     "Normalize linear 1-in-1-out segments once after loading GFA");
    convert.add_flag("--segment-no-id", parse.convert_args.segment_no_id,
                     "Store BGFA segments without explicit segment id (id assigned by read order)");
    convert.add_flag("--no-sosr", parse.convert_args.no_sosr,
                     "Disable SOSR/SN storage in BGFA segments");
    convert.add_flag("--sosr-compact", parse.convert_args.sosr_compact,
                     "Use compact SOSR+SN packing in BGFA segments");

    // merge
    auto &merge = *app.add_subcommand("merge", "merge some gfas/bgfas to one");
    merge.add_option("-g,--graph", parse.merge_args.graph_path,
                     "Input GFA/BGFA file")
        ->expected(-1)
        ->required()
        ->check(CLI::ExistingFile);
    merge.add_option("-o,--output", parse.merge_args.output_path,
                     "Output file name")
        ->default_str("merged.bgfa");
    merge.add_option("-t,--type", parse.merge_args.merge_type,
                     "Merge type (0: strict, 1: exactly)")
        ->default_val(0)
        ->check(CLI::Range(0, 1));
    merge.add_option("-l,--merge-limited", parse.merge_args.merge_limited,
                     "Merge similarity threshold for exact merge (1-100)")
        ->default_val(90)
        ->check(CLI::Range(1, 100));
    // new options
    merge.add_option("-k,--alignment-kernel",
                     parse.merge_args.alignment_kernel,
                     "Alignment kernel (wfa, ksw)")
        ->default_str("wfa")
        ->transform(CLI::CheckedTransformer(merge_alig_map, CLI::ignore_case));
    merge.add_option("-d,--details-type", parse.merge_args.details_type,
                     "Auxiliary method (0:none, 1:distance, 2:minimizer)")
        ->default_val(0)
        ->check(CLI::Range(0, 2));

    // index
    auto &index = *app.add_subcommand("index", "Create index for GFA/BGFA");
    auto index_group = index.add_option_group("Input file");
    index_group->add_option("-g,--gfa", parse.index_args.input_path,
                            "Input GFA file")
        ->check(CLI::ExistingFile);
    index_group->add_option("-b,--bgfa", parse.index_args.input_path,
                            "Input BGFA file")
        ->check(CLI::ExistingFile);
    index_group->require_option(1); // 必须指定其中一个选项
    index.add_option("-m, --mode", parse.index_args.seed_mode,
                     "Seed mode (kmer, minimizer)")
        ->default_str("minimizer")
        ->check(CLI::IsMember({"kmer", "minimizer"}));
    index.add_option("--minimizer-mode", parse.index_args.iterate_mode,
                     "Indexing mode: (single, single-reverse, double)\n\t\tsingle: Iterate over all nodes in the graph only, without considering any base complementary sequences\n\t\tsingle-reverse: Iterate over all nodes in the graph and each k-mer will prefer the kmer with the smallest hash in its and its base complementary sequences\n\t\tdouble: Iterate over all nodes and their base complements")
        ->default_str("double")
        ->check(CLI::IsMember({0, 1, 2}))
        ->transform(CLI::CheckedTransformer(index_itrmode_map,
                                            CLI::ignore_case));
    index.add_option("-w,--window", parse.index_args.window_size,
                     "Minimizer Window size")
        ->default_val(5)
        ->check(CLI::Range(1, 100));
    index.add_option("-k,--kmer-length", parse.index_args.minimizer_length,
                     "k-mer length")
        ->default_val(15)
        ->check(CLI::Range(1, 32));
    index.add_option("-o,--output", parse.index_args.output_path,
                     "Output BMIN file")
        ->default_str("auto");

    // view
    auto &view = *app.add_subcommand("view", "Show information");
    view.add_option("-g,--bgfa", parse.view_args.graph_path, "bGFA file")
        ->check(CLI::ExistingFile);
    view.add_option("-m,--bmin", parse.view_args.minimizer_path, "bMIN file")
        ->check(CLI::ExistingFile);

    // stats
    auto &stats = *app.add_subcommand("stats", "Show graph statistics");
    stats.add_option("-g,--bgfa", parse.view_args.graph_path, "GFA/BGFA file")
        ->required()
        ->check(CLI::ExistingFile);

    CLI11_PARSE(app, argc, argv);

    if (convert)
    {
        time = ConvertProcess(parse);
    }
    else if (index)
    {
        // if (index.count("--gfa"))
        //     parse.index_args.is_bgfa = false;
        // else if (index.count("--bgfa"))
        //     parse.index_args.is_bgfa = true;
        time = IndexProcess(parse);
    }
    else if (view)
    {
        time = ViewProcess(parse);
    }
    else if (merge)
    {
        time = MergeProcess(parse);
    }
    else if (stats)
    {
        time = StatsProcess(parse);
    }
    else
    {
        std::cerr << "Error: No valid command provided." << std::endl;
        return 0;
    }
    if (parse.time)
    {
        std::cout << "Total time: " << time << " seconds." << std::endl;
    }

    return 0;
}