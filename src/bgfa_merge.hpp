/*
 * @file \bGFA\src\bgfa_merge.hpp
 * @brief
 * @details
 * @author grampart
 * @version
 * @copyright: Copyright(c) 2025 by grampart, All Rights Reserved.
 * @license: MIT License
 */
#ifndef BGFA_MERGE_HPP_
#define BGFA_MERGE_HPP_

#include "bgfa.hpp"
#include "bgfa_graph.hpp"
#include "bgfa_seed.hpp"
#include "bgfa_alignment.hpp"

/// Placeholder size, used for graph ID shifting in merging
#define MERGE_PLACEHOLDER 4
/// Placeholder size, used for renumbering segments in exactly merging
#define EXACT_MERGE_PLACEHOLDER 12

/// Dynamic calculation for clearing the EXACT_MERGE_PLACEHOLDER bit field mask
#define EXACT_MERGE_CLEAR_MASK (uword_t(-1) ^ (((uword_t(1) << EXACT_MERGE_PLACEHOLDER) - 1) << (WORD_BIT - EXACT_MERGE_PLACEHOLDER)))

#define SCREEN_MODE_DIRECT 0    ///< All nodes are considered for merging
#define SCREEN_MODE_MINIMIZER 2 ///< Nodes are screened based on minimizers
#define SCREEN_MODE_DISTANCE 4  ///< Nodes are screened based on distance
#define SCREEN_MODE_PATH 8      ///< Nodes are screened based on paths
#define MERGE_MODE_STRICT 0     ///< Merge nodes if they are consistent
#define MERGE_MODE_EXACT 1      ///< Merge nodes based on alignment

/// New ID mapping: graph_id | old_id -> (new_first_id, new_last_id)
typedef std::unordered_map<uword_t, std::pair<uword_t, uword_t>> NewIdMap;

// TODO: 目前仅设计了增量合并，是否需要增加归并合并？

/**
 * @brief Update links in the merged graph based on segment set
 * @param {GFA} &merged_graph - The merged GFA graph
 * @param {vector<GFA>} &graphs - List of GFA graphs to merge
 * @param {NewIdMap} &new_id_list - Mapping of graph_id | old_id to new ID ranges
 */
void updateLink(GFA &merged_graph,
                const std::vector<GFA> &graphs,
                const NewIdMap &new_id_list)
{
    // 0 is the reference graph
    uword_t graph_id = 1;

    for (const auto &gfa : graphs)
    {
        const auto &all_links = gfa.getLinksSet().getAllLinks();
        for (uword_t from_id = 0; from_id < all_links.size(); from_id++)
        {
            const auto &data = all_links[from_id];
            for (const auto &to_value : data)
            {

                uword_t graph_from_id =
                    ((graph_id << (WORD_BIT - MERGE_PLACEHOLDER)) | from_id);

                uword_t graph_to_id =
                    ((graph_id << (WORD_BIT - MERGE_PLACEHOLDER)) |
                     (to_value >> LINK_DIR_BIT));

                auto from_it = new_id_list.find(graph_from_id);
                auto to_it = new_id_list.find(graph_to_id);

                if (from_it != new_id_list.end() && to_it != new_id_list.end())
                {
                    // From last id -> To first id
                    merged_graph.addLink(
                        from_it->second.second,
                        (to_it->second.first << LINK_DIR_BIT | (to_value & 0x03)));
                }
            }
        }

        graph_id++;
    }
    return;
}

/**
 * @brief Update segment set with segments from a GFA graph
 * @param {NodeMap} &segment_set - Mapping of segment sequences to new IDs
 * @param {GFA} &gfa - GFA graph to extract segments from
 * @param {uword_t} graph_id - Unique ID for the graph
 */
void updateSegmentSet(NodeMap &segment_set, const GFA &gfa, uword_t graph_id)
{
    for (const auto &segment : gfa.getAllSegments())
    {
        std::string segment_seq = segment.getSequenceAsString();
        segment_set[segment_seq] =
            graph_id << (WORD_BIT - MERGE_PLACEHOLDER) | segment.getId();
    }
}

/**
 * @brief Ensure a sequence exists in merged_graph and segment_set, return its index
 * If the sequence is already present in segment_set, return the mapped index.
 * Otherwise, create a new Segment in merged_graph, record it in segment_set and return the new index.
 */
uword_t ensureSegmentInMerged(GFA &merged_graph, NodeMap &segment_set, const std::string &seq)
{
    auto it = segment_set.find(seq);
    if (it != segment_set.end())
    {
        return it->second;
    }
    Segment tmp(0, seq);
    uword_t idx = merged_graph.addSegment(tmp);
    segment_set[seq] = idx;
    return idx;
}

uword_t getSimilarScore(MinimizerSet &seg_ms, MinimizerSet &ref_ms)
{
    MinimizerSetVect seg_ms_vect = getMinimizerVect(seg_ms);
    MinimizerSetVect ref_ms_vect = getMinimizerVect(ref_ms);

    size_t m = seg_ms_vect.size();
    size_t n = ref_ms_vect.size();
    if (m == 0 || n == 0)
        return 0;

    // TODO: 硬编码，后续换成参数传入
    const int MATCH = 1;       // hash match +1
    const int ORDER_BONUS = 1; // additional +1 if previous match preserves order
    const int MISMATCH = -1;   // mismatch -1
    const int GAP = -1;        // gap penalty -1

    // Use rolling 1D DP (two rows) for Smith-Waterman local alignment to save memory.
    std::vector<int> prev(n + 1, 0);
    std::vector<int> curr(n + 1, 0);
    int max_score = 0;

    for (size_t i = 1; i <= m; ++i)
    {
        curr[0] = 0;
        for (size_t j = 1; j <= n; ++j)
        {
            int diag_score;
            // compare hashes (second element of pair)
            if (seg_ms_vect[i - 1].second == ref_ms_vect[j - 1].second)
            {
                diag_score = prev[j - 1] + MATCH;
                // if previous items exist and offsets are increasing in both sequences,
                // treat this as preserving order and give an extra bonus
                if (i > 1 && j > 1 &&
                    seg_ms_vect[i - 2].first < seg_ms_vect[i - 1].first &&
                    ref_ms_vect[j - 2].first < ref_ms_vect[j - 1].first)
                {
                    diag_score += ORDER_BONUS;
                }
            }
            else
            {
                diag_score = prev[j - 1] + MISMATCH;
            }

            int up_score = prev[j] + GAP;
            int left_score = curr[j - 1] + GAP;

            int cell = std::max(0, std::max(diag_score, std::max(up_score, left_score)));
            curr[j] = cell;
            if (cell > max_score)
                max_score = cell;
        }

        // move current row to previous, and reset current
        std::swap(prev, curr);
        std::fill(curr.begin(), curr.end(), 0);
    }

    return (uword_t)max_score;
}

uword_t getMergeTargetbyMinimizers(GFA &merged_graph,
                                   MinimizerSet &ref_ms, MinimizerSet &seg_ms)
{
    uword_t target_seg_id = uword_t(-1);
    // id -> (current_score, last_seg_offset)
    std::unordered_map<uword_t, uword_t> seg_minimizer_hit;

    // hit count for each minimizer in the segment
    for (const auto &seg_minimizer : seg_ms.getMinimizerTable())
    {
        uint64_t minimizer_hash = seg_minimizer.first;

        // Check if this minimizer exists in the reference minimizer set
        auto ref_it = ref_ms.getMinimizerTable().find(minimizer_hash);

        // Calculate hits
        if (ref_it != ref_ms.getMinimizerTable().end())
        {
            for (const auto &pos : ref_it->second.getPositions())
            {
                uword_t ref_seg_id = getPosID(pos);
                uword_t ref_seg_offset = getPosOffset(pos);

                auto score_it = seg_minimizer_hit.find(ref_seg_id);
                if (score_it != seg_minimizer_hit.end())
                {
                    seg_minimizer_hit[ref_seg_id] += 1;
                }
                else
                {
                    seg_minimizer_hit[ref_seg_id] = 1;
                }
            }
        }
    }

    if (seg_minimizer_hit.empty())
    {
        return target_seg_id;
    }

    // hit order
    // TODO: 阈值不应该是硬编码，后续需要调整
    uword_t min_hit_count = seg_ms.getMinimizerTable().size() * 0.7;
    uword_t max_score = 0;
    for (const auto &hit_pair : seg_minimizer_hit)
    {
        if (hit_pair.second < min_hit_count)
        {
            continue;
        }

        target_seg_id = hit_pair.first;
        Segment target_seg = merged_graph.getSegmentByIndex(hit_pair.first);
        MinimizerSet target_seg_ms = MinimizerSet(
            target_seg,
            ref_ms.getLength(),
            ref_ms.getWindowSize(),
            ref_ms.getMode());
        if (max_score < getSimilarScore(seg_ms, target_seg_ms))
        {
            max_score = getSimilarScore(seg_ms, target_seg_ms);
            target_seg_id = hit_pair.first;
        }
    }
    return target_seg_id;
}

/**
 * @brief Merge a segment with a target sequence using exact comparison
 * @param {GFA} &merged_graph - The merged GFA graph
 * @param {NewIdMap} &new_id_list - Mapping of graph_id | old_id to new ID
 * @param {uword_t} g_id - Unique ID for the graph
 * @param {Segment} &segment - Segment to be merged
 * @param {string} &target_seq - Target sequence to compare against
 * @param {uword_t} target_id - Target segment ID in merged_graph
 * @param {bool} &is_merged - Output flag indicating if the segment was merged
 * @param {uword_t} merge_limited - Merge threshold percentage
 * @param {uint8_t} alignment_kernel - Alignment kernel to use (0: KSW, 1: WFA)
 * @return {*}
 */
void mergeExactCompare(GFA &merged_graph, NewIdMap &new_id_list, uword_t g_id,
                       const Segment &segment, const std::string &target_seq,
                       const uword_t target_id, bool &is_merged,
                       uword_t merge_limited, uint8_t alignment_kernel)
{
    const std::string segment_seq = segment.getSequenceAsString();
    const uword_t segment_id = segment.getId();

    // Skip if length difference exceeds threshold
    double length_diff = std::abs((int)segment_seq.size() - (int)target_seq.size());
    length_diff /= std::max(segment_seq.size(), target_seq.size());
    if (length_diff > (100 - merge_limited) / 100.0)
        return;

#ifdef INFO_DEBUG
    std::cerr << "[mergeExactCompare] segment_id=" << segment_id
              << " target_id=" << target_id
              << " seg_len=" << segment_seq.size()
              << " target_len=" << target_seq.size()
              << " seg_seq=" << segment_seq
              << " target_seq=" << target_seq
              << " length_diff=" << length_diff << "\n";
#endif

    std::string cigar_str = alignmentSequence(segment_seq, target_seq);
    if (cigar_str.empty() || cigar_str == "0M")
        return;

    uint64_t m_number = 0;
    uint64_t cigar_number = 0;
    uint64_t total_ops = 0;
    {
        size_t p = 0, n = cigar_str.size();
        const char *s1 = segment_seq.data();
        const char *s2 = target_seq.data();
        size_t pos1_tmp = 0, pos2_tmp = 0;
        while (p < n)
        {
            uint64_t num = 0;
            while (p < n && std::isdigit(static_cast<unsigned char>(cigar_str[p])))
            {
                num = num * 10 + (cigar_str[p] - '0');
                ++p;
            }
            if (p >= n)
                break; // malformed
            char op = cigar_str[p++];
            ++total_ops;
            cigar_number += num;

            if (op == 'M')
            {

                m_number += num;
                pos1_tmp += num;
                pos2_tmp += num;
            }
            else if (op == 'I')
            {
                pos2_tmp += num;
            }
            else if (op == 'D')
            {
                pos1_tmp += num;
            }
            else if (op == 'X')
            {
                pos1_tmp += num;
                pos2_tmp += num;
            }
        }
    }

    if (cigar_number == 0)
        return;

#ifdef INFO_DEBUG
    std::cerr << "\t\t    cigar=" << cigar_str
              << " cigar_len=" << cigar_number
              << " m_number=" << m_number
              << " total_ops=" << total_ops
              << " merge_limited=" << merge_limited << "\n";
#endif

    if (m_number * 100 < merge_limited * cigar_number)
        return;
    if (total_ops > (static_cast<uint64_t>(1) << EXACT_MERGE_PLACEHOLDER) - 1)
        return;

    // Default: not merged
    is_merged = false;

    uword_t first_new_target = -1, last_new_target = -1;
    uword_t first_new_segment = -1, last_new_segment = -1;

    size_t pos1 = 0, pos2 = 0;
    bool is_first_seg1 = true, is_first_seg2 = true;
    uword_t seq1_cur_id = uword_t(-1), seq2_cur_id = uword_t(-1);

    {
        size_t p = 0, n = cigar_str.size();
        uword_t op_idx = 0;
        const char *s1 = segment_seq.data();
        const char *s2 = target_seq.data();

        while (p < n)
        {
            uint64_t num = 0;
            while (p < n && std::isdigit(static_cast<unsigned char>(cigar_str[p])))
            {
                num = num * 10 + (cigar_str[p] - '0');
                ++p;
            }
            if (p >= n)
                break;
            char op = cigar_str[p++];

            uword_t id_tag = static_cast<uword_t>(op_idx) << (WORD_BIT - EXACT_MERGE_PLACEHOLDER);
            ++op_idx;

            if (op == 'M')
            {
                pos1 += num;
                pos2 += num;
                const size_t start1 = pos1 - num;
                const size_t start2 = pos2 - num;

                // create (or reuse) merged segments for the matched region
                if (num > 0)
                {
                    std::string piece1 = segment_seq.substr(start1, num);
                    std::string piece2 = target_seq.substr(start2, num);

                    if (piece1 == piece2)
                    {
                        Segment tmp(0, piece1);
                        uword_t new_idx = merged_graph.addSegment(tmp);
                        if (is_first_seg1)
                        {
                            first_new_segment = new_idx;
                            is_first_seg1 = false;
                        }
                        last_new_segment = new_idx;

                        if (is_first_seg2)
                        {
                            first_new_target = new_idx;
                            is_first_seg2 = false;
                        }
                        last_new_target = new_idx;

                        if (seq1_cur_id != uword_t(-1))
                        {
                            merged_graph.addLink(seq1_cur_id, (new_idx << LINK_DIR_BIT) | 0x03);
                        }
                        seq1_cur_id = new_idx;

                        if (seq2_cur_id != uword_t(-1))
                        {
                            merged_graph.addLink(seq2_cur_id, (new_idx << LINK_DIR_BIT) | 0x03);
                        }
                        seq2_cur_id = new_idx;
                    }
                    else
                    {
                        // different sequences in M -> create two nodes
                        Segment tmp1(0, piece1);
                        Segment tmp2(0, piece2);
                        uword_t new_idx1 = merged_graph.addSegment(tmp1);
                        uword_t new_idx2 = merged_graph.addSegment(tmp2);

                        if (is_first_seg1)
                        {
                            first_new_segment = new_idx1;
                            is_first_seg1 = false;
                        }
                        last_new_segment = new_idx1;

                        if (is_first_seg2)
                        {
                            first_new_target = new_idx2;
                            is_first_seg2 = false;
                        }
                        last_new_target = new_idx2;

                        if (seq1_cur_id != uword_t(-1))
                        {
                            merged_graph.addLink(seq1_cur_id, (new_idx1 << LINK_DIR_BIT) | 0x03);
                        }
                        seq1_cur_id = new_idx1;

                        if (seq2_cur_id != uword_t(-1))
                        {
                            merged_graph.addLink(seq2_cur_id, (new_idx2 << LINK_DIR_BIT) | 0x03);
                        }
                        seq2_cur_id = new_idx2;
                    }
                }
            }
            else if (op == 'X')
            {
                pos1 += num;
                pos2 += num;
                const size_t start1 = pos1 - num;
                const size_t start2 = pos2 - num;

                // Create separate nodes for mismatched region on both sides
                if (num > 0)
                {
                    std::string piece1 = segment_seq.substr(start1, num);
                    std::string piece2 = target_seq.substr(start2, num);

                    Segment tmp1(0, piece1);
                    Segment tmp2(0, piece2);
                    uword_t new_idx1 = merged_graph.addSegment(tmp1);
                    uword_t new_idx2 = merged_graph.addSegment(tmp2);

                    if (is_first_seg1)
                    {
                        first_new_segment = new_idx1;
                        is_first_seg1 = false;
                    }
                    last_new_segment = new_idx1;
                    if (seq1_cur_id != uword_t(-1))
                    {
                        merged_graph.addLink(seq1_cur_id, (new_idx1 << LINK_DIR_BIT) | 0x03);
                    }
                    seq1_cur_id = new_idx1;

                    if (is_first_seg2)
                    {
                        first_new_target = new_idx2;
                        is_first_seg2 = false;
                    }
                    last_new_target = new_idx2;
                    if (seq2_cur_id != uword_t(-1))
                    {
                        merged_graph.addLink(seq2_cur_id, (new_idx2 << LINK_DIR_BIT) | 0x03);
                    }
                    seq2_cur_id = new_idx2;
                }
            }
            else if (op == 'I')
            {
                pos2 += num;
                const size_t start2 = pos2 - num;
                if (num > 0)
                {
                    std::string piece2 = target_seq.substr(start2, num);
                    Segment tmp2(0, piece2);
                    uword_t new_idx2 = merged_graph.addSegment(tmp2);
                    if (is_first_seg2)
                    {
                        first_new_target = new_idx2;
                        is_first_seg2 = false;
                    }
                    last_new_target = new_idx2;
                    if (seq2_cur_id != uword_t(-1))
                    {
                        merged_graph.addLink(seq2_cur_id, (new_idx2 << LINK_DIR_BIT) | 0x03);
                    }
                    seq2_cur_id = new_idx2;
                }
            }
            else if (op == 'D')
            {
                pos1 += num;
                const size_t start1 = pos1 - num;
                if (num > 0)
                {
                    std::string piece1 = segment_seq.substr(start1, num);
                    Segment tmp1(0, piece1);
                    uword_t new_idx1 = merged_graph.addSegment(tmp1);
                    if (is_first_seg1)
                    {
                        first_new_segment = new_idx1;
                        is_first_seg1 = false;
                    }
                    last_new_segment = new_idx1;
                    if (seq1_cur_id != uword_t(-1))
                    {
                        merged_graph.addLink(seq1_cur_id, (new_idx1 << LINK_DIR_BIT) | 0x03);
                    }
                    seq1_cur_id = new_idx1;
                }
            }
            else
            {
                continue;
            }
        }
    }

    // Only consider "merged" if new nodes are actually created on the source side, and write to the mapping
    if (first_new_segment != uword_t(-1) && last_new_segment != uword_t(-1))
    {
        is_merged = true;

        // Redirect links for the Target segment（仅当目标侧也产生了首尾时）
        if (first_new_target != uword_t(-1) && last_new_target != uword_t(-1))
        {
            // Redirect all incoming edges X -> target to X -> first_new_target
            merged_graph.redirectLinksTo(target_id, first_new_target);

            // Move all outgoing edges target -> Y to last_new_target -> Y
            merged_graph.redirectLinksFrom(target_id, last_new_target);

            // Remove any remaining edges of target (if any) and mark segment deleted
            merged_graph.removeLinksOf(target_id);
            merged_graph.markSegmentDeleted(target_id);
        }

        // Set Information in new_id_list
        new_id_list[(g_id << (WORD_BIT - MERGE_PLACEHOLDER)) | segment_id] =
            std::make_pair(first_new_segment, last_new_segment);
    }
}

/**
 * @brief Merge multiple GFA graphs strictly by matching node sequences
 * @details Merges nodes with identical sequences across graphs. Only exact
 *          matches are considered; no sequence similarity is used. The
 *          resulting graph may be larger than with other methods if inputs
 *          have many similar but non-identical nodes.
 * @param {GFA} &merged_graph - Reference to the GFA
 * @param {vector<GFA>} graphs - Vector of GFA graphs to merge
 * @return {void}
 */
void mergeGFAStrict(GFA &merged_graph, const std::vector<GFA> &graphs)
{
    uword_t graph_id = 0;
    NodeMap segment_set;  ///< segment sequence -> new segment ID
    NewIdMap new_id_list; ///< old segment ID -> new segment ID

    updateSegmentSet(segment_set, merged_graph, 0);
    graph_id++;

    for (const auto &gfa : graphs)
    {
        for (const auto &segment : gfa.getAllSegments())
        {
            std::string segment_seq = segment.getSequenceAsString();

            if (segment_set.find(segment_seq) == segment_set.end())
            {
                segment_set[segment_seq] = merged_graph.addSegment(segment);
            }
            else
            {
                new_id_list[(graph_id << (WORD_BIT - MERGE_PLACEHOLDER)) |
                            segment.getId()] =
                    std::make_pair(segment_set[segment_seq],
                                   segment_set[segment_seq]);
            }
        }
        graph_id++;
    }

    updateLink(merged_graph, graphs, new_id_list);
}

void mergeGFAExact(GFA &merged_graph,
                   const std::vector<GFA> &graphs,
                   uword_t merge_limited,
                   uint8_t alignment_kernel)
{
    uword_t graph_id = 0;
    NodeMap segment_set;  ///< segment sequence -> new segment ID (only for not merged segments)
    NewIdMap new_id_list; ///< old segment ID -> new segment ID

    if (merge_limited <= 0 || merge_limited > 100)
    {
        throw std::logic_error(
            "[Merge GFA Exact] merge_limited should be between 1 and 100.");
    }

    updateSegmentSet(segment_set, merged_graph, 0);
    graph_id++;

    for (const auto &gfa : graphs)
    {
        for (const auto &segment : gfa.getAllSegments())
        {
            std::string segment_seq = segment.getSequenceAsString();
            uword_t segment_id = segment.getId();

            if (segment_set.find(segment_seq) == segment_set.end())
            {
                bool is_merged = false;

                // std::string segment_seq = segment.getSequenceAsString();
                // uword_t segment_id = segment.getId();

                for (const auto &segment_in_set : segment_set)
                {
                    if (is_merged)
                        break;
                    mergeExactCompare(merged_graph, new_id_list, graph_id,
                                      segment, segment_in_set.first,
                                      segment_in_set.second, is_merged,
                                      merge_limited, alignment_kernel);
                }
                if (is_merged)
                {
#ifdef INFO_DEBUG
                    std::cerr << "[mergeGFAExact] Merged segment_id="
                              << segment_id << " seq=" << segment_seq << "\n";
#endif
                    continue;
                }
                else
                {
                    uword_t new_id = merged_graph.addSegment(segment);

                    // Record new segment ID in new_id_list and segment_set
                    new_id_list[(graph_id << (WORD_BIT - MERGE_PLACEHOLDER)) | segment_id] =
                        std::make_pair(new_id, new_id);
                    segment_set[segment_seq] = new_id;
#ifdef INFO_DEBUG
                    std::cerr
                        << "[mergeGFAExact] New segment_id="
                        << segment_id << " seq=" << segment_seq << "\n";
#endif
                }
            }
            else
            {
                new_id_list[(graph_id << (WORD_BIT - MERGE_PLACEHOLDER)) |
                            segment.getId()] =
                    std::make_pair(segment_set[segment_seq],
                                   segment_set[segment_seq]);
            }
        }
        graph_id++;
    }

    updateLink(merged_graph, graphs, new_id_list);
}

/**
 * @brief Merge multiple GFA graphs strictly by distance
 * @details This function merges graphs based on segment distances.
 *          Segments from different graphs that have similar distances
 *          (within a defined threshold) are considered for merging.
 * @note This function may be enabled, as it appears that distance
 *       optimizations in strict merges are not faster than a single
 *       find in a set?
 * @param {GFA} &merged_graph
 * @param {vector<GFA>} &graphs
 * @return {*}
 */
void mergeGFAStrictByDistance(GFA &merged_graph, const std::vector<GFA> &graphs)
{
    uword_t graph_id = 0;
    NodeMap segment_set;
    NewIdMap new_id_list;

    updateSegmentSet(segment_set, merged_graph, 0);
    graph_id++;

    for (const auto &gfa : graphs)
    {
        for (const auto &segment : gfa.getAllSegments())
        {
            uword_t seg_dis = segment.getDis();
            std::vector<Segment> candidates; // Const &
            // TODO: 此处阈值需要设置代码进行传入(DISTANCE_THRESHOLD)，而不是硬编码
            candidates = merged_graph.findSegmentsByDistance(seg_dis, 10);
            bool is_found = false;
            for (const auto &cand_seg : candidates)
            {
                if (cand_seg.getSequenceAsString() ==
                    segment.getSequenceAsString())
                {
                    std::string segment_seq = segment.getSequenceAsString();
                    new_id_list[(graph_id << (WORD_BIT - MERGE_PLACEHOLDER)) |
                                segment.getId()] =
                        std::make_pair(segment_set[segment_seq],
                                       segment_set[segment_seq]);

                    is_found = true;
                    break;
                }
            }
            if (!is_found)
            {
                segment_set[segment.getSequenceAsString()] = merged_graph.addSegment(segment);
            }
        }
        graph_id++;
    }

    updateLink(merged_graph, graphs, new_id_list);
}

void mergeGFAExactByDistance(GFA &merged_graph,
                             const std::vector<GFA> &graphs,
                             uword_t merge_limited,
                             uint8_t alignment_kernel)
{
    uword_t graph_id = 0;
    NodeMap segment_set;
    NewIdMap new_id_list;

    if (merge_limited <= 0 || merge_limited > 100)
    {
        throw std::logic_error(
            "[Merge GFA Exact] merge_limited should be between 1 and 100.");
    }

    updateSegmentSet(segment_set, merged_graph, 0);
    graph_id++;
    for (const auto &gfa : graphs)
    {
        for (const auto &segment : gfa.getAllSegments())
        {
            std::string segment_seq = segment.getSequenceAsString();
            uword_t segment_id = segment.getId();

            if ((segment_set.find(segment_seq) == segment_set.end()) ||
                segment_set.empty())
            {
                bool is_merged = false;

                // Find candidates by distance
                uword_t seg_dis = segment.getDis();
                std::vector<Segment> candidates;
                // TODO: 此处阈值需要设置代码进行传入(DISTANCE_THRESHOLD)，而不是硬编码
                candidates = merged_graph.findSegmentsByDistance(seg_dis, 10);

                for (const auto &cand_seg : candidates)
                {
                    if (is_merged)
                        break;
                    mergeExactCompare(merged_graph, new_id_list, graph_id,
                                      segment, cand_seg.getSequenceAsString(),
                                      cand_seg.getId(), is_merged,
                                      merge_limited, alignment_kernel);
                }

                if (is_merged)
                {
#ifdef INFO_DEBUG
                    std::cerr << "[mergeGFAExact] Merged segment_id="
                              << segment_id << " seq=" << segment_seq << "\n";
#endif
                    continue;
                }
                else
                {
                    uword_t new_id = merged_graph.addSegment(segment);

                    // Record new segment ID in new_id_list and segment_set
                    new_id_list[(graph_id << (WORD_BIT - MERGE_PLACEHOLDER)) | segment_id] =
                        std::make_pair(new_id, new_id);
                    segment_set[segment_seq] = new_id;
#ifdef INFO_DEBUG
                    std::cerr << "[mergeGFAExact] New segment_id="
                              << segment_id << " seq=" << segment_seq << "\n";
#endif
                }
            }
            else
            {
                new_id_list[(graph_id << (WORD_BIT - MERGE_PLACEHOLDER)) |
                            segment.getId()] =
                    std::make_pair(segment_set[segment_seq],
                                   segment_set[segment_seq]);
            }
        }
        graph_id++;
    }

    updateLink(merged_graph, graphs, new_id_list);
}

/**
 * @brief Merge multiple GFA graphs exactly by minimizers
 * @details This function generates default empty MinimizerSet instances for
 *          each graph and then calls mergeGFAStrictByMinimize() that utilizes
 *          minimizers for merging
 * @param {GFA} &merged_graph
 * @param {vector<GFA>} &graphs
 * @return {*}
 */
void mergeGFAExactByMinimize(GFA &merged_graph,
                             const std::vector<GFA> &graphs,
                             uword_t merge_limited,
                             uint8_t alignment_kernel)
{
    // TODO: minimizer的长度和窗口大小需要作为参数传入，目前仅是硬编码
    MinimizerSet ref_ms = MinimizerSet(merged_graph, 15, 20);
    uword_t graph_id = 0;
    NodeMap segment_set;
    NewIdMap new_id_list;

    if (merge_limited <= 0 || merge_limited > 100)
    {
        throw std::logic_error(
            "[Merge GFA Exact] merge_limited should be between 1 and 100.");
    }

    updateSegmentSet(segment_set, merged_graph, 0);
    graph_id++;

    for (const auto &graph : graphs)
    {
        // MinimizerSet temp_ms = MinimizerSet(graph, 15, 20);
        for (const auto &segment : graph.getAllSegments())
        {
            std::string segment_seq = segment.getSequenceAsString();
            uword_t segment_id = segment.getId();

            if (segment_set.find(segment_seq) == segment_set.end())
            {
                bool is_merged = false;
                // Minimizer of a segment
                MinimizerSet seg_ms = MinimizerSet(segment,
                                                   ref_ms.getLength(),
                                                   ref_ms.getWindowSize(),
                                                   ref_ms.getMode());

                uword_t target_id = getMergeTargetbyMinimizers(merged_graph,
                                                               ref_ms,
                                                               seg_ms);

                // Guard: if no valid target found, add as new segment
                if (target_id == uword_t(-1))
                {
                    segment_set[segment_seq] = merged_graph.addSegment(segment);
                }
                else
                {
                    std::string target_seq =
                        merged_graph.getSegmentByIndex(target_id).getSequenceAsString();
                    mergeExactCompare(merged_graph, new_id_list, graph_id,
                                      segment, target_seq, target_id, is_merged,
                                      merge_limited, alignment_kernel);

                    if (!is_merged)
                    {
                        segment_set[segment_seq] = merged_graph.addSegment(segment);
                    }
                }
            }
            else
            {
                new_id_list[(graph_id << (WORD_BIT - MERGE_PLACEHOLDER)) |
                            segment.getId()] =
                    std::make_pair(segment_set[segment_seq],
                                   segment_set[segment_seq]);
            }
        }
        graph_id++;
    }

    updateLink(merged_graph, graphs, new_id_list);
}

/**
 * @brief Merge multiple GFA graphs fuzzily
 * @param {vector<GFA>} graphs - Vector of GFA graphs to merge
 * @return {GFA} Merged GFA graph
 */
GFA mergeGraphs(const GFA &ref_graph,
                const std::vector<GFA> &graphs,
                uword_t merge_args,
                uword_t merge_limited = 90,
                uint8_t alignment_kernel = 1)
{
    GFA merged_graph = ref_graph;

    if (graphs.empty())
    {
        throw std::invalid_argument("[Merge Graphs] Input graph list is empty");
    }

    // TODO: merge_args 需要包含更多信息，例如minimizer长度、窗口大小等参数
    switch (merge_args)
    {
    case SCREEN_MODE_DIRECT | MERGE_MODE_STRICT:
        mergeGFAStrict(merged_graph, graphs);
        break;
    case SCREEN_MODE_MINIMIZER | MERGE_MODE_STRICT:
        std::cout << "Due to performance considerations,"
                  << "strict merging using the minimizer has been abandoned."
                  << "We will use exact merging as an alternative.";
        mergeGFAExactByMinimize(
            merged_graph, graphs, merge_limited, alignment_kernel);
        break;
    case SCREEN_MODE_DISTANCE | MERGE_MODE_STRICT:
        mergeGFAStrictByDistance(merged_graph, graphs);
        break;
    case SCREEN_MODE_PATH | MERGE_MODE_STRICT:
        // mergeGFAStrictByPath(merged_graph, graphs);
        break;
    case SCREEN_MODE_DIRECT | MERGE_MODE_EXACT:
        mergeGFAExact(merged_graph, graphs, merge_limited, alignment_kernel);
        break;
    case SCREEN_MODE_MINIMIZER | MERGE_MODE_EXACT:
        mergeGFAExactByMinimize(
            merged_graph, graphs, merge_limited, alignment_kernel);
        break;
    case SCREEN_MODE_DISTANCE | MERGE_MODE_EXACT:
        mergeGFAExactByDistance(merged_graph, graphs, merge_limited, alignment_kernel);
        break;
    case SCREEN_MODE_PATH | MERGE_MODE_EXACT:
        // mergeGFAExactByPath(merged_graph, graphs);
        break;

    default:
        break;
    }

    merged_graph.compactAndRenumber();

    return merged_graph;
}

#endif // MERGE_HPP_