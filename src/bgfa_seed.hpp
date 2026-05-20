/*
 * @file \bGFA\src\bgfa_seed.hpp
 * @brief
 * @details
 * @author grampart
 * @version
 * @copyright: Copyright(c) 2025 by grampart, All Rights Reserved.
 * @license: MIT License
 */
#ifndef BGFA_SEED_HPP_
#define BGFA_SEED_HPP_

#include "bgfa.hpp"
#include "bgfa_graph.hpp"

typedef uint64_t uhash_t;
typedef uint64_t upos_t;
typedef uint64_t useq_t;

typedef std::vector<std::pair<uword_t, uhash_t>> MinimizerSetVect;

#define MINIMIZER_MAX_LENGTH 32 ///< Maximum length of minimizer

#define MINIMIZER_ID_SHITF 32            ///< Bits for segment ID
#define MINIMIZER_ID_MASK 0xFFFFFFFF     ///< Mask for segment ID
#define MINIMIZER_OFFSET_SHITF 1         ///< Bits for offset
#define MINIMIZER_OFFSET_MSAK 0x7FFFFFFF ///< Mask for offset
#define MINIMIZER_DIR_SHITF 0            ///< Bits for direction
#define MINIMIZER_DIR_MASK 0x01          ///< Mask for direction

#define MINIMIZER_MODE_SINGLE 0         ///< Single strand minimizer mode
#define MINIMIZER_MODE_SINGLE_REVERSE 1 ///< Single strand with reverse complement minimizer mode
#define MINIMIZER_MODE_DOUBLE 2         ///< Double strand minimizer mode

/**
 * @brief 64-bit hash function
 * @details Implements a simple 64-bit hash function for generating hash values
 *          from keys. The original source of this code is
 *          https://github.com/lh3/minigraph/blob/master/sketch.c
 * @param {uint64_t} key
 * @param {uint64_t} mask
 * @return {uint64_t} Hashed value
 */
uint64_t hash64(uint64_t key, uint64_t mask)
{
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

// 似乎需要一个新的类表示Minimizer?

uword_t getPosID(upos_t pos)
{
    return (uword_t)(pos >> MINIMIZER_ID_SHITF & MINIMIZER_ID_MASK);
}

uword_t getPosOffset(upos_t pos)
{
    return (uword_t)((pos >> MINIMIZER_OFFSET_SHITF) & MINIMIZER_OFFSET_MSAK);
}

uword_t getPosDir(upos_t pos)
{
    return (uword_t)((pos >> MINIMIZER_DIR_SHITF) & MINIMIZER_DIR_MASK);
}

class Minimizer
{
private:
    useq_t sequence_;               ///< Minimizer hash value
    std::vector<upos_t> positions_; ///< Positions of the minimizer

public:
    /**
     * @brief Default constructor initializes an empty Minimizer
     */
    Minimizer() : sequence_(0) {}

    Minimizer(useq_t seq) : sequence_(seq) {}

    Minimizer(useq_t seq, const std::vector<upos_t> &pos)
        : sequence_(seq), positions_(pos) {}

    void addPosition(upos_t pos) { positions_.push_back(pos); }

    void setSequence(useq_t seq)
    {
        if (positions_.size() > 0 && sequence_ != seq)
        {
            throw std::runtime_error("[Minimizer Set Sequence] Error: Cannot change sequence of existing minimizer.");
        }
        sequence_ = seq;
    }

    const useq_t &getSequence() const { return sequence_; }

    const std::vector<upos_t> &getPositions() const { return positions_; }

    friend std::ostream &operator<<(std::ostream &os, Minimizer &ms)
    {
        os << "Minimizer Sequence: " << ms.sequence_ << ", Positions: ";
        for (const auto &pos : ms.positions_)
        {
            os << "\t id: " << getPosID(pos)
               << ", offset: " << getPosOffset(pos)
               << ", dir:: " << getPosDir(pos) << ";\n";
        }
        return os;
    }
};

/**
 * @brief Class for managing minimizers in a graph
 * @details This class provides methods to extract and store minimizers from graph segments.
 */
class MinimizerSet
{
private:
    uword_t length_;      ///< Length of the minimizer
    uword_t window_size_; ///< Window size for minimizer selection
    uword_t mode_;        ///< Minimizer mode (single, single-reverse, double)

    /**
     * @brief Minimizer table (hash map)
     * @details key->hash; value->(sequence(64), [offset(32+31+1)])
     *          offset: 32-bits id, 31-bits offset,1 bit direction
     */
    std::unordered_map<uhash_t, Minimizer> hash_minimizer_map_;

public:
    /**
     * @brief Default constructor initializes with default minimizer length and window size
     * @details Default length is 15, window size is 20 and mode is single strand
     */
    MinimizerSet() : length_(15), window_size_(20), mode_(MINIMIZER_MODE_SINGLE) {}

    /**
     * @brief Constructor with specified minimizer length and window size
     * @param {int} length_ - Length of the minimizer
     * @param {int} window_size_ - Window size for minimizer selection
     */
    explicit MinimizerSet(int length, int window_size) : length_(length), window_size_(window_size), mode_(MINIMIZER_MODE_SINGLE)
    {
        if (length > 32)
        {
            throw std::runtime_error("[MinimizerSet Constructor] Error: length exceeds 32.");
        }
    }

    /**
     * @brief Constructor with specified minimizer length, window size, and mode
     * @param {int} length_ - Length of the minimizer
     * @param {int} window_size_ - Window size for minimizer selection
     * @param {uint8_t} mode - Minimizer mode
     */
    explicit MinimizerSet(uword_t length, uword_t window_size, uword_t mode)
        : length_(length), window_size_(window_size), mode_(mode)
    {
        if (length > MINIMIZER_MAX_LENGTH)
        {
            throw std::runtime_error("[MinimizerSet Constructor] Error: length exceeds" + std::to_string(MINIMIZER_MAX_LENGTH) + ".");
        }
    }

    /**
     * @brief Constructor that extracts minimizers from a GFA graph
     * @param {GFA} gfa - GFA graph object
     * @param {int} length_ - Length of the minimizer
     * @param {int} window_size_ - Window size for minimizer selection
     */
    explicit MinimizerSet(GFA gfa, uword_t length, uword_t window_size)
        : MinimizerSet(length, window_size)
    {
        getMinimizerFromGFA(gfa);
    }

    /**
     * @brief Constructor that extracts minimizers from a GFA graph with specified mode
     * @param {GFA} gfa - GFA graph object
     * @param {int} length_ - Length of the minimizer
     * @param {int} window_size_ - Window size for minimizer selection
     */
    explicit MinimizerSet(
        GFA gfa, uword_t length, uword_t window_size, uword_t mode)
        : MinimizerSet(length, window_size, mode)
    {
        getMinimizerFromGFA(gfa);
    }

    /**
     * @brief Constructor that extracts minimizers from a Segment
     * @param {Segment} segment - Segment object
     * @param {int} length_ - Length of the minimizer
     * @param {int} window_size_ - Window size for minimizer selection
     */
    explicit MinimizerSet(Segment segment, uword_t length, uword_t window_size)
        : MinimizerSet(length, window_size)
    {
        getMinimizerFromSegment(segment);
    }

    /**
     * @brief Constructor that extracts minimizers from a Segment with specified mode
     * @param {Segment} segment - Segment object
     * @param {int} length_ - Length of the minimizer
     * @param {int} window_size_ - Window size for minimizer selection
     */
    explicit MinimizerSet(
        Segment segment, uword_t length, uword_t window_size, uword_t mode)
        : MinimizerSet(length, window_size, mode)
    {
        getMinimizerFromSegment(segment);
    }

    /**
     * @brief Constructor that loads minimizers from a bMIN file
     * @param {string} &filename - Path to the bMIN file
     */
    explicit MinimizerSet(const std::string &filename)
    {
        if (filename.size() < 5 || filename.substr(filename.size() - 5) != ".bmin")
        {
            throw std::runtime_error(
                "[MinimizerSet Constructor] Filename is not a .bmin file: " +
                filename);
        }

        getMinimizerFrombmin(filename);
    }

    /**
     * @brief Get the number of minimizers stored
     * @return {uword_t} Number of minimizers
     */
    uword_t getMinimizerNum() const { return hash_minimizer_map_.size(); }

    /**
     * @brief Get the minimizer length
     * @return {uword_t} Minimizer length
     */
    uword_t getLength() const { return length_; }

    /**
     * @brief Get the minimizer window size
     * @return {uword_t} Minimizer window size
     */
    uword_t getWindowSize() const { return window_size_; }

    /**
     * @brief Get the minimizer mode
     * @return {uword_t} Minimizer mode
     */
    uword_t getMode() const { return mode_; }

    /**
     * @brief Get the minimizer table
     * @return {const map<uint64_t, pair<uint64_t, vector<uint64_t>>>&} Minimizer table
     */
    const std::unordered_map<uhash_t, Minimizer> &getMinimizerTable() const
    {
        return hash_minimizer_map_;
    }

    /**
     * @brief Set the minimizer length, it can only be set when the minimizer table is empty
     * @param {uword_t} length - New minimizer length
     */
    void setLength(uword_t length)
    {
        if (hash_minimizer_map_.size() > 0)
        {
            throw std::runtime_error("[Set Minimizer Length] Error: Cannot set length when minimizer table is not empty.");
        }
        if (length > MINIMIZER_MAX_LENGTH)
        {
            throw std::runtime_error("[Set Minimizer Length] Error: length exceeds" + std::to_string(MINIMIZER_MAX_LENGTH) + ".");
        }
        length_ = length;
    }

    /**
     * @brief Set the minimizer window size, it can only be set when the minimizer table is empty
     * @param {uword_t} window_size - New minimizer window size
     */
    void setWindowSize(uword_t window_size)
    {
        if (hash_minimizer_map_.size() > 0)
        {
            throw std::runtime_error("[Set Minimizer Window Size] Error: Cannot set window size when minimizer table is not empty.");
        }
        window_size_ = window_size;
    }

    /**
     * @brief Set the minimizer mode, it can only be set when the minimizer table is empty
     * @param {uword_t} mode - New minimizer mode
     */
    void setMode(uword_t mode)
    {
        if (hash_minimizer_map_.size() > 0)
        {
            throw std::runtime_error("[Set Minimizer Mode] Error: Cannot set mode when minimizer table is not empty.");
        }
        mode_ = mode;
    }

    /**
     * @brief Output stream operator for GFA
     */
    friend std::ostream &operator<<(std::ostream &os, MinimizerSet &ms)
    {
        os << "Minimizer Length: " << std::dec << ms.getLength() << std::endl;
        os << "Minimizer Window Size: " << std::dec << ms.getWindowSize() << std::endl;

        for (const auto &it : ms.getMinimizerTable())
        {
            std::cout << "\tHash Value: 0x" << std::hex << it.first;
            std::cout << "\tSequence: " << ms.bin2sequence(it.second.getSequence(), ms.getLength()) << std::endl;
            for (int i = 0; i < it.second.getPositions().size(); i++)
            {
                std::cout << "\t\tid:0x" << std::hex
                          << std::setw(8) << std::setfill('0')
                          << getPosID(it.second.getPositions()[i]);
                std::cout << "\toffset:0x" << std::hex
                          << std::setw(8) << std::setfill('0')
                          << getPosOffset(it.second.getPositions()[i]);
                std::cout << "\tdir:" << std::hex
                          << ((getPosDir((it.second.getPositions()[i])) == 0) ? '+' : '-')
                          << std::endl;
            }
        }
    }

    /**
     * @brief Print minimizer data in binary format
     * @details Outputs the minimizer length, window size, and each minimizer entry in hexadecimal format.
     *          It often used for debugging or exporting minimizer data.
     */
    void print2bin()
    {
        std::cout << "0x" << std::hex << std::setw(8) << std::setfill('0') << length_ << " ";
        std::cout << "0x" << std::hex << std::setw(8) << std::setfill('0') << window_size_ << std::endl;
        for (auto it : hash_minimizer_map_)
        {
            for (int i = 0; i < it.second.getPositions().size(); i++)
            {
                std::cout << "0x" << std::hex << std::setw(16) << std::setfill('0') << it.first << " ";
                std::cout << "0x" << std::hex << std::setw(16) << std::setfill('0') << it.second.getSequence() << " ";
                std::cout << "0x" << std::hex << std::setw(16) << std::setfill('0') << it.second.getPositions()[i] << std::endl;
            }
            std::cout << std::endl;
        }
    }

    void writeSeed2bmin(const std::string &filename)
    {
        std::ofstream output_file(filename, std::ios::binary);
        if (!output_file.is_open())
        {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }

        output_file.write(
            reinterpret_cast<const char *>(&length_), sizeof(uint32_t));
        output_file.write(
            reinterpret_cast<const char *>(&window_size_), sizeof(uint32_t));
        for (auto it : hash_minimizer_map_)
        {
            for (int i = 0; i < it.second.getPositions().size(); i++)
            {
                output_file.write(reinterpret_cast<const char *>(&(it.first)), sizeof(uint64_t));
                output_file.write(reinterpret_cast<const char *>(&(it.second.getSequence())), sizeof(uint64_t));
                output_file.write(reinterpret_cast<const char *>(&(it.second.getPositions()[i])), sizeof(uint64_t));
            }
        }

        output_file.close();
    }

private:
    /**
     * @brief Get minimizers from a graph
     * @param {GFA} gfa
     */
    void getMinimizerFromGFA(GFA gfa)
    {
        for (auto segment : gfa.getAllSegments())
        {
            // std::cout << "node:" << node[0] << std::endl;
            getMinimizerFromSegment(segment);
        }
    }

    /**
     * @brief Load minimizers from a bMIN file
     * @param {string} &filename - Path to the bMIN file
     */
    void getMinimizerFrombmin(const std::string &filename)
    {
        std::ifstream input_file(filename, std::ios::binary);
        if (!input_file.is_open())
        {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }

        // 读取length_和window_size_
        uword_t read_length, read_window;
        input_file.read(
            reinterpret_cast<char *>(&read_length), sizeof(uword_t));
        input_file.read(
            reinterpret_cast<char *>(&read_window), sizeof(uword_t));
        length_ = read_length;
        window_size_ = read_window;

        if (length_ > MINIMIZER_MAX_LENGTH)
        {
            throw std::runtime_error(
                "[Load From bmin File] Error: length_ exceeds 32.");
            // std::cerr << "Error: length_ exceeds 32." << std::endl;
            input_file.close();
            return;
        }

        uhash_t hash;
        useq_t seq;
        upos_t pos;

        while (input_file.read(
            reinterpret_cast<char *>(&hash), sizeof(uhash_t)))
        {
            if (!input_file.read(
                    reinterpret_cast<char *>(&seq), sizeof(useq_t)))
            {
                throw std::runtime_error(
                    "[Load From bmin File] Error reading sequence.");
            }
            if (!input_file.read(
                    reinterpret_cast<char *>(&pos), sizeof(upos_t)))
            {
                throw std::runtime_error(
                    "[Load From bmin File] Error reading position.");
            }

            auto &entry = hash_minimizer_map_[hash];
            if (entry.getPositions().empty())
            {
                entry.setSequence(seq); // 首次遇到该哈希，设置序列
                // entry.first = seq;
            }
            else if (entry.getSequence() != seq)
            {
                throw std::runtime_error(
                    "[Load From bmin File] Conflict: Hash 0x" +
                    std::to_string(hash) + " has inconsistent sequences (0x" +
                    std::to_string(entry.getSequence()) + " vs 0x" +
                    std::to_string(seq) + ")");
                // std::cerr << "Conflict: Hash 0x" << std::hex << hash
                //           << " has inconsistent sequences (0x" << entry.first
                //           << " vs 0x" << seq << ")" << std::endl;
                // continue; // 跳过不一致的条目
            }
            entry.addPosition(pos);
            // entry.second.push_back(pos);
        }

        input_file.close();
        std::cout << "Loaded " << hash_minimizer_map_.size() << " minimizers from " << filename << std::endl;
    }

    /**
     * @brief Convert binary encoded sequence to string
     * @param {uint64_t} bin_strm - Binary encoded sequence
     * @param {uint16_t} l_str - Length of the sequence
     * @return {string} Decoded sequence string
     */
    std::string bin2sequence(useq_t bin_strm, uword_t l_str)
    {
        uword_t bit_number = 2 * l_str;
        std::string sequence;
        for (size_t i = 0; i < l_str; ++i)
        {
            uint8_t baseBinary = (bin_strm >> (bit_number - 2 * i - 2)) & 0x03;
            sequence += GlobalVariant::bin2char[baseBinary];
        }
        return sequence;
    }

    /**
     * @brief Get minimizers from a segment based on the specified mode
     * @param {Segment} segment
     */
    void getMinimizerFromSegment(Segment segment)
    {
        switch (mode_)
        {
        case MINIMIZER_MODE_SINGLE:
            getMinimizerOfNode_SingleMode(segment);
            break;
        case MINIMIZER_MODE_SINGLE_REVERSE:
            getMinimizerOfNode_SingleReverseMode(segment);
            break;
        case MINIMIZER_MODE_DOUBLE:
            getMinimizerOfNode_DoubleMode(segment);
            break;
        default:
            std::cerr << "Error: Invalid minimizer mode." << std::endl;
            break;
        }
    }

    /**
     * @brief Get minimizers of a node in single strand mode
     * @details The core concept of this function is based on Li Heng's minimap2 code on GitHub (https://github.com/lh3/minimap2)
     * @param {Segment} node
     * @param {bool} is_reverse
     * @return {*}
     */
    void getMinimizerOfNode_SingleMode(Segment node, bool is_reverse = false)
    {
        int k = length_;
        int w = window_size_;
        uint64_t shift1 = 2 * (k - 1);
        uint64_t mask = (1ULL << 2 * k) - 1;
        uint64_t kmer = 0;
        int i, j, l, buf_pos, min_pos;
        // int kmer_span = 0;

        // index信息, [0]:hash, [1]:node_id(32),offset(31),direction(1), [2]:sequence
        uint64_t min[3] = {UINT64_MAX, UINT64_MAX, 0};
        uint64_t buf[256][3];

        int l_str = node.getLength();
        uint8_t baseinElement = WORD_BIT / 2;

        for (i = l = buf_pos = min_pos = 0; i < l_str; ++i)
        {
            int index_seq_element = i / baseinElement;
            int index_base = i % baseinElement;

            assert(index_seq_element < node.getBinarySequence().size());

            uint8_t c = (node.getBinarySequence()[index_seq_element] >> (WORD_BIT - 2 * index_base - 2)) & 0x03;
            uint64_t hashinfo[3] = {UINT64_MAX, UINT64_MAX, 0};

            if (is_reverse)
                kmer = (kmer >> 2) | (3ULL ^ c) << shift1;
            else
                kmer = (kmer << 2 | c) & mask;
            ++l;

            if (l >= k)
            {
                hashinfo[0] = hash64(kmer, mask);
                if (is_reverse)
                    hashinfo[1] = (uint64_t)node.getId() << 32 | (l_str + k - 2 - (uint32_t)i) << 1 | 1;
                else
                    hashinfo[1] = (uint64_t)node.getId() << 32 | (uint32_t)i << 1 | 0;
                hashinfo[2] = kmer;
            }

            buf[buf_pos][0] = hashinfo[0];
            buf[buf_pos][1] = hashinfo[1];
            buf[buf_pos][2] = hashinfo[2];

            if (l == w + k - 1 && min[0] != UINT64_MAX)
            {
                for (j = buf_pos + 1; j < w; ++j)
                {
                    if (min[0] == buf[j][0] && buf[j][1] != min[1])
                    {
                        hash_minimizer_map_[buf[j][0]].setSequence(buf[j][2]);
                        hash_minimizer_map_[buf[j][0]].addPosition(buf[j][1]);
                        // hash_minimizer_map_[buf[j][0]].first = buf[j][2];
                        // hash_minimizer_map_[buf[j][0]].second.emplace_back(buf[j][1]);
                        // info_table[buf[j][0]][1]++;
                    }
                }
                for (j = 0; j < buf_pos; ++j)
                {
                    if (min[0] == buf[j][0] && buf[j][1] != min[1])
                    {
                        hash_minimizer_map_[buf[j][0]].setSequence(buf[j][2]);
                        hash_minimizer_map_[buf[j][0]].addPosition(buf[j][1]);

                        // hash_minimizer_map_[buf[j][0]].first = buf[j][2];
                        // hash_minimizer_map_[buf[j][0]].second.emplace_back(buf[j][1]);
                        // hash_minimizer_map_[buf[j][0]][1]++;
                    }
                }
            }
            if (hashinfo[0] <= min[0])
            {
                if (l >= w + k && min[0] != UINT64_MAX)
                {
                    hash_minimizer_map_[min[0]].setSequence(min[2]);
                    hash_minimizer_map_[min[0]].addPosition(min[1]);

                    // hash_minimizer_map_[min[0]].first = min[2];
                    // hash_minimizer_map_[min[0]].second.emplace_back(min[1]);
                    // hash_minimizer_map_[min[0]][1]++;
                }
                min[0] = hashinfo[0];
                min[1] = hashinfo[1];
                min[2] = hashinfo[2];
                min_pos = buf_pos;
            }
            else if (buf_pos == min_pos)
            {
                if (l >= w + k - 1 && min[0] != UINT64_MAX)
                {
                    hash_minimizer_map_[min[0]].setSequence(min[2]);
                    hash_minimizer_map_[min[0]].addPosition(min[1]);

                    // hash_minimizer_map_[min[0]].first = min[2];
                    // hash_minimizer_map_[min[0]].second.emplace_back(min[1]);
                    // hash_minimizer_map_[min[0]][1]++;
                }
                min[0] = UINT64_MAX;
                for (j = buf_pos + 1; j < w; ++j)
                {
                    if (min[0] >= buf[j][0])
                    {
                        min[0] = buf[j][0];
                        min[1] = buf[j][1];
                        min[2] = buf[j][2];
                        min_pos = j;
                    }
                }
                for (j = 0; j <= buf_pos; ++j)
                {
                    if (min[0] >= buf[j][0])
                    {
                        min[0] = buf[j][0];
                        min[1] = buf[j][1];
                        min[2] = buf[j][2];
                        min_pos = j;
                    }
                }

                if (l >= w + k - 1 && min[0] != UINT64_MAX)
                {
                    for (j = buf_pos + 1; j < w; ++j)
                    {
                        if (min[0] == buf[j][0] && min[1] != buf[j][1])
                        {
                            hash_minimizer_map_[buf[j][0]].setSequence(buf[j][2]);
                            hash_minimizer_map_[buf[j][0]].addPosition(buf[j][1]);

                            // hash_minimizer_map_[buf[j][0]].first = buf[j][2];
                            // hash_minimizer_map_[buf[j][0]].second.emplace_back(buf[j][1]);
                            // hash_minimizer_map_[buf[j][0]][1]++;
                        }
                    }
                    for (j = 0; j <= buf_pos; ++j)
                    {
                        if (min[0] == buf[j][0] && min[1] != buf[j][1])
                        {
                            hash_minimizer_map_[buf[j][0]].setSequence(buf[j][2]);
                            hash_minimizer_map_[buf[j][0]].addPosition(buf[j][1]);

                            // hash_minimizer_map_[buf[j][0]].first = buf[j][2];
                            // hash_minimizer_map_[buf[j][0]].second.emplace_back(buf[j][1]);
                            // hash_minimizer_map_[buf[j][0]][1]++;
                        }
                    }
                }
            }
            if (++buf_pos == w)
                buf_pos = 0;
        }

        if (min[0] != UINT64_MAX)
        {
            hash_minimizer_map_[min[0]].setSequence(min[2]);
            hash_minimizer_map_[min[0]].addPosition(min[1]);

            // hash_minimizer_map_[min[0]].first = min[2];
            // hash_minimizer_map_[min[0]].second.emplace_back(min[1]);
            // hash_minimizer_map_[min[0]][1]++;
        }
    }

    /**
     * @brief Get minimizers of a node in single strand reverse mode
     * @details The core concept of this function is based on Li Heng's minimap2 code on GitHub (https://github.com/lh3/minimap2)
     * @param {Segment} node
     * @return {*}
     */
    void getMinimizerOfNode_SingleReverseMode(Segment node)
    {
        int k = length_;
        int w = window_size_;
        uint64_t shift1 = 2 * (k - 1);
        uint64_t mask = (1ULL << 2 * k) - 1;
        uint64_t kmer[2] = {0, 0};
        int i, j, l, buf_pos, min_pos;
        // int kmer_span = 0;

        // index信息, [0]:hash, [1]:node_id(32),offset(31),direction(1), [2]:sequence
        uint64_t min[3] = {UINT64_MAX, UINT64_MAX, 0};
        uint64_t buf[256][3];

        int l_str = node.getLength();
        uint8_t baseinElement = WORD_BIT / 2;

        for (i = l = buf_pos = min_pos = 0; i < l_str; ++i)
        {
            int index_seq_element = i / baseinElement;
            int index_base = i % baseinElement;

            assert(index_seq_element < node.getBinarySequence().size());

            uint8_t c = (node.getBinarySequence()[index_seq_element] >> (WORD_BIT - 2 * index_base - 2)) & 0x03;
            uint64_t hashinfo[3] = {UINT64_MAX, UINT64_MAX, 0};
            uint8_t kmer_dir;

            kmer[0] = (kmer[0] << 2 | c) & mask;
            kmer[1] = (kmer[1] >> 2) | (3ULL ^ c) << shift1;
            if (kmer[0] == kmer[1])
                continue;
            kmer_dir = kmer[0] < kmer[1] ? 0 : 1;
            ++l;

            if (l >= k)
            {
                hashinfo[0] = hash64(kmer[kmer_dir], mask);
                hashinfo[1] = (uint64_t)node.getId() << 32 | (uint32_t)i << 1 | kmer_dir;
                hashinfo[2] = kmer[kmer_dir];
            }

            buf[buf_pos][0] = hashinfo[0];
            buf[buf_pos][1] = hashinfo[1];
            buf[buf_pos][2] = hashinfo[2];

            if (l == w + k - 1 && min[0] != UINT64_MAX)
            {
                for (j = buf_pos + 1; j < w; ++j)
                {
                    if (min[0] == buf[j][0] && buf[j][1] != min[1])
                    {
                        hash_minimizer_map_[buf[j][0]].setSequence(buf[j][2]);
                        hash_minimizer_map_[buf[j][0]].addPosition(buf[j][1]);

                        // hash_minimizer_map_[buf[j][0]].first = buf[j][2];
                        // hash_minimizer_map_[buf[j][0]].second.emplace_back(buf[j][1]);
                        // info_table[buf[j][0]][1]++;
                    }
                }
                for (j = 0; j < buf_pos; ++j)
                {
                    if (min[0] == buf[j][0] && buf[j][1] != min[1])
                    {
                        hash_minimizer_map_[buf[j][0]].setSequence(buf[j][2]);
                        hash_minimizer_map_[buf[j][0]].addPosition(buf[j][1]);

                        // hash_minimizer_map_[buf[j][0]].first = buf[j][2];
                        // hash_minimizer_map_[buf[j][0]].second.emplace_back(buf[j][1]);
                        // hash_minimizer_map_[buf[j][0]][1]++;
                    }
                }
            }
            if (hashinfo[0] <= min[0])
            {
                if (l >= w + k && min[0] != UINT64_MAX)
                {
                    hash_minimizer_map_[min[0]].setSequence(min[2]);
                    hash_minimizer_map_[min[0]].addPosition(min[1]);

                    // hash_minimizer_map_[min[0]].first = min[2];
                    // hash_minimizer_map_[min[0]].second.emplace_back(min[1]);
                    // // hash_minimizer_map_[min[0]][1]++;
                }
                min[0] = hashinfo[0];
                min[1] = hashinfo[1];
                min[2] = hashinfo[2];
                min_pos = buf_pos;
            }
            else if (buf_pos == min_pos)
            {
                if (l >= w + k - 1 && min[0] != UINT64_MAX)
                {
                    hash_minimizer_map_[min[0]].setSequence(min[2]);
                    hash_minimizer_map_[min[0]].addPosition(min[1]);

                    // hash_minimizer_map_[min[0]].first = min[2];
                    // hash_minimizer_map_[min[0]].second.emplace_back(min[1]);
                    // hash_minimizer_map_[min[0]][1]++;
                }
                min[0] = UINT64_MAX;
                for (j = buf_pos + 1; j < w; ++j)
                {
                    if (min[0] >= buf[j][0])
                    {
                        min[0] = buf[j][0];
                        min[1] = buf[j][1];
                        min[2] = buf[j][2];
                        min_pos = j;
                    }
                }
                for (j = 0; j <= buf_pos; ++j)
                {
                    if (min[0] >= buf[j][0])
                    {
                        min[0] = buf[j][0];
                        min[1] = buf[j][1];
                        min[2] = buf[j][2];
                        min_pos = j;
                    }
                }

                if (l >= w + k - 1 && min[0] != UINT64_MAX)
                {
                    for (j = buf_pos + 1; j < w; ++j)
                    {
                        if (min[0] == buf[j][0] && min[1] != buf[j][1])
                        {
                            hash_minimizer_map_[buf[j][0]].setSequence(buf[j][2]);
                            hash_minimizer_map_[buf[j][0]].addPosition(buf[j][1]);

                            // hash_minimizer_map_[buf[j][0]].first = buf[j][2];
                            // hash_minimizer_map_[buf[j][0]].second.emplace_back(buf[j][1]);
                            // hash_minimizer_map_[buf[j][0]][1]++;
                        }
                    }
                    for (j = 0; j <= buf_pos; ++j)
                    {
                        if (min[0] == buf[j][0] && min[1] != buf[j][1])
                        {
                            hash_minimizer_map_[buf[j][0]].setSequence(buf[j][2]);
                            hash_minimizer_map_[buf[j][0]].addPosition(buf[j][1]);
                            // hash_minimizer_map_[buf[j][0]].first = buf[j][2];
                            // hash_minimizer_map_[buf[j][0]].second.emplace_back(buf[j][1]);
                            // hash_minimizer_map_[buf[j][0]][1]++;
                        }
                    }
                }
            }
            if (++buf_pos == w)
                buf_pos = 0;
        }

        if (min[0] != UINT64_MAX)
        {
            hash_minimizer_map_[min[0]].setSequence(min[2]);
            hash_minimizer_map_[min[0]].addPosition(min[1]);

            // hash_minimizer_map_[min[0]].first = min[2];
            // hash_minimizer_map_[min[0]].second.emplace_back(min[1]);
            // hash_minimizer_map_[min[0]][1]++;
        }
    }

    /**
     * @brief Get minimizers of a node in double strand mode
     * @details It calls both single strand and single strand reverse mode functions
     * @param {Segment} node
     * @return {*}
     */
    void getMinimizerOfNode_DoubleMode(Segment node)
    {
        getMinimizerOfNode_SingleMode(node, false);
        getMinimizerOfNode_SingleMode(node, true);
    }
};

void testMinimizerSet(MinimizerSet ms)
{
    std::cout << "Minimizer Total number: " << std::dec << ms.getMinimizerTable().size() << std::endl;
    std::cout << "txt" << std::endl;
    std::cout << ms;
    std::cout << "bin" << std::endl;
    ms.print2bin();
}

MinimizerSetVect getMinimizerVect(MinimizerSet &ms)
{
    MinimizerSetVect minimizer_vect;
    for (const auto &it : ms.getMinimizerTable())
    {
        uhash_t hash = it.first;
        for (auto pos : it.second.getPositions())
        {
            uword_t node_offset = getPosOffset(pos);
            minimizer_vect.emplace_back(std::make_pair(node_offset, hash));
        }
    }

    std::sort(minimizer_vect.begin(), minimizer_vect.end());

    return minimizer_vect;
}

#endif // SEED_HPP_