/**
 * @file /bGFA/src/graph.hpp
 * @brief Graph data structures and utilities for bGFA format
 * @details This file defines the core data structures and functions for
 *          handling graphs in the bGFA format, including segments and links. *          It includes four main classes: Segment, LinksSet, Paths and Graph.
 * @author grampart
 * @version V1.0
 * @copyright Copyright(c) 2025 by grampart, All Rights Reserved.
 * @license MIT License
 */

#ifndef BGFA_GRAPH_HPP_
#define BGFA_GRAPH_HPP_

#define _MEMORY_ALIGN true

#include "bgfa.hpp"
#include <cmath>
#include <numeric>
#include <map>
#include <random>
#include <limits>

/// Placeholder size, used for representing link direction
#define LINK_DIR_BIT 2
/// Placeholder size, used for flag/length in compact segme nt representation
#define SEGMENT_COMPACT_PLACEHOLDER 6
/// Placeholder size, used for offset/rank in compact segment representation
#define SOSR_OFFSET_SHIFT 6
/// Mask for rank in sosr
#define SR_MASK 0x3F
/// Maximum increment for path encoding
#define PATH_MAX_INC 0xFF
/// Shift for increment encoding in path
#define PATH_INC_SHIFT 8

/// pre_info bit0: path mode (0 direct, 1 increment)
#define PRE_INFO_PATHMODE_INC 0x01
/// pre_info bit1: segment id storage (0 store id, 1 no id in segment record)
#define PRE_INFO_SEGMENT_NO_ID 0x02
/// pre_info bit2: SOSR stored (0 no, 1 yes)
#define PRE_INFO_HAS_SOSR 0x04
/// pre_info bit3: SOSR storage mode (0 compressed, 1 uncompressed)
#define PRE_INFO_SOSR_UNCOMPRESSED 0x08

#if WORD_BIT == 16
#define PATH_INC_SIZE 2
#elif WORD_BIT == 32
#define PATH_INC_SIZE 4
#elif WORD_BIT == 64
#define PATH_INC_SIZE 8
#endif

namespace GlobalVariant
{
    std::array<uint8_t, 256> char2bin = []
    {
        std::array<uint8_t, 256> arr{};
        arr.fill(255); // 填充无效值
        arr['A'] = arr['a'] = 0;
        arr['C'] = arr['c'] = 1;
        arr['G'] = arr['g'] = 2;
        arr['T'] = arr['t'] = 3;
        return arr;
    }();

    // 二进制到字符的映射
    std::array<char, 4> bin2char = {'A', 'C', 'G', 'T'};

    // 方向字符到二进制的映射
    std::array<uint8_t, 256> dir2bin = []
    {
        std::array<uint8_t, 256> arr{};
        arr.fill(255);
        arr['+'] = 1;
        arr['-'] = 0;
        return arr;
    }();

    // 二进制到方向字符
    std::array<char, 2> bin2dir = {'-', '+'};

    // CIGAR 操作符到二进制
    std::array<uint8_t, 256> cigarop2bin = []
    {
        std::array<uint8_t, 256> arr{};
        arr.fill(255);
        arr['M'] = 0;
        arr['I'] = 1;
        arr['D'] = 2;
        arr['S'] = 3;
        arr['N'] = 4;
        arr['H'] = 5;
        arr['P'] = 6;
        arr['='] = 7;
        arr['X'] = 8;
        arr['B'] = 9;
        return arr;
    }();

    // 二进制到 CIGAR 操作符
    std::array<char, 10> bin2cigarop =
        {'M', 'I', 'D', 'S', 'N', 'H', 'P', '=', 'X', 'B'};
}

typedef std::unordered_map<std::string, uword_t> NodeMap;

enum PathMode
{
    PATHMODE_DIRECT = 0,
    PATHMODE_INCREMENT = 1,
    PATHMODE_INCREMENTHYBRID = 2
};

/**
 * @brief Segment class representing a graph segment
 * @details The equality comparison for segments does not follow the same rules
 *          as the greater-than and less-than comparisons. Two segments are
 *          considered equal if their sequences are identical, regardless of
 *          their IDs or other attributes. A segment is considered less than
 *          another if its dis is smaller, and then by id.
 */
class Segment
{
private:
    uword_t id_;                    ///< Segment ID
    uword_t dis_;                   ///< Segment distance, not stored
    uword_t so_;                    ///< Segment offset (SO)
    uword_t sr_;                    ///< Segment rank (SR)
    uword_t sn_id_;                 ///< Segment name id (SN)
    bool has_sosr_;                 ///< Whether SOSR/SN are stored
    bool sosr_compact_;             ///< SOSR/SN compact storage mode
    uword_t length_;                ///< Segment length(SL:i), not stored
    std::vector<uword_t> sequence_; ///< Binary encoded sequence

    uword_t store_size_; ///< Size in storage (in uword_t units)

public:
    /**
     * @brief Default constructor initializes an empty segment
     */
    Segment() : so_(0), sr_(0), sn_id_(0), has_sosr_(true), sosr_compact_(false), length_(0), store_size_(0) {}

    /**
     * @brief Constructor with id, dis, length, and binary sequence
     * @param {uword_t} id - Segment ID
     * @param {uword_t} sosr - Segment offset and rank
     * @param {uword_t} length - Segment length
     * @param {vector<uword_t>} &sequence - Binary encoded sequence
     */
    explicit Segment(const uword_t id, uword_t so, uword_t sr, uword_t length,
                     const std::vector<uword_t> &sequence)
        : id_(id), so_(so), sr_(sr), sn_id_(0), has_sosr_(true), sosr_compact_(false),
          length_(length), sequence_(sequence)
    {
        setStoreSize();
    }

    /**
     * @brief Constructor with id, length, and binary sequence (default: 0)
     * @param {uword_t} id - Segment ID
     * @param {uword_t} length - Segment length
     * @param {vector<uword_t>} &sequence - Binary encoded sequence
     */
    explicit Segment(const uword_t id, uword_t length,
                     const std::vector<uword_t> &sequence)
        : id_(id), so_(0), sr_(0), sn_id_(0), has_sosr_(true), sosr_compact_(false),
          length_(length), sequence_(sequence), dis_(0)
    {
        setStoreSize();
    }

    /**
     * @brief Constructor with id and string sequence (dis default: 0)
     * @param {uword_t} id - Segment ID
     * @param {string} &sequence_str - Segment sequence string
     */
    explicit Segment(const uword_t id, const std::string &sequence_str)
        : id_(id), so_(0), sr_(0), sn_id_(0), has_sosr_(true), sosr_compact_(false),
          length_(sequence_str.length()), dis_(0)
    {
        stringToBinary(sequence_str);
        setStoreSize();
    }

    explicit Segment(const uword_t id, const std::string &sequence_str, const uword_t offset, const uword_t rank)
        : id_(id), so_(0), sr_(0), sn_id_(0), has_sosr_(true), sosr_compact_(false),
          length_(sequence_str.length()), dis_(0)
    {
        setSosr(offset, rank);
        stringToBinary(sequence_str);
        setStoreSize();
    }

    /**
     * @brief Get the Id of the Segment
     * @return {uword_t} Segment ID
     */
    uword_t getId() const { return id_; }

    /**
     * @brief Get the Dis of the Segment
     * @return {uword_t} Segment distance
     */
    uword_t getDis() const { return dis_; }

    /**
     * @brief Get the SOSR of the Segment
     * @return {uword_t} Segment sosr
     */
    uword_t getSosr() const
    {
        if (!has_sosr_)
            return 0;
        if (sosr_compact_)
        {
            const uword_t sr = (sr_ & 0x07);
            const uword_t so = (so_ & 0x003FFFFF);
            const uword_t sn = (sn_id_ & 0x7F);
            return (sn << 25) | (so << 3) | sr;
        }
        const uword_t sr = (sr_ & 0x7F);
        const uword_t so = (so_ & 0x01FFFFFF);
        return (so << 7) | sr;
    }

    /**
     * @brief Get the Offset of the Segment
     * @return {uword_t} Segment offset
     */
    uword_t getOffset() const { return has_sosr_ ? so_ : 0; }

    /**
     * @brief Get the Rank of the Segment
     * @return {uword_t} Segment rank
     */
    uword_t getRank() const { return has_sosr_ ? sr_ : 0; }

    uword_t getSnId() const { return sn_id_; }

    bool hasSosr() const { return has_sosr_; }

    bool isSosrCompact() const { return sosr_compact_; }

    /**
     * @brief Get the Length of the Segment
     * @return {uword_t} Segment length
     */
    uword_t getLength() const { return length_; }

    /**
     * @brief Get the storage size of the Segment in uword_t units
     * @return {uword_t} Storage size
     */
    uword_t getStoreSize()
    {
        setStoreSize();
        return store_size_;
    }

    uword_t getStoreSize(bool store_id) const
    {
        uword_t sz = 0;
        if (store_id)
            sz += 1; // id
        if (has_sosr_)
            sz += sosr_compact_ ? 1 : 2; // sosr (+sn if uncompressed)

        if (store_id)
        {
            if (length_ > (WORD_BIT - 6) / 2)
                sz += (1 + static_cast<uword_t>(sequence_.size()));
            else
                sz += 1;
            return sz;
        }

        const uword_t compact1_cap = (WORD_BIT - 6) / 2;
        const uword_t compact1_len_max = std::min<uword_t>(compact1_cap, 0x0F);
        const uword_t compact2_cap = (2 * WORD_BIT - 8) / 2;
        const uword_t compact2_len_max = std::min<uword_t>(compact2_cap, 0x3F);

        if (length_ > compact2_len_max)
            sz += (1 + static_cast<uword_t>(sequence_.size()));
        else if (length_ <= compact1_len_max)
            sz += 1;
        else
            sz += 2;

        return sz;
    }

    /**
     * @brief Get the Sequence of the Segment as a string
     * @return {string} Segment sequence
     */
    std::string getSequenceAsString() const
    {
        std::string sequence;
        for (size_t i = 0; i < length_; ++i)
        {
            sequence += getBaseAt(i);
        }
        return sequence;
    }

    /**
     * @brief Get the binary encoded Sequence of the Segment
     * @return {const vector<uword_t>&} Binary encoded sequence
     */
    const std::vector<uword_t> &getBinarySequence() const { return sequence_; }

    /**
     * @brief Get the base character at a specific position in the Segment
     * @param {size_t} position - Position in the sequence
     * @return {char} Base character at the specified position
     */
    char getBaseAt(size_t position) const
    {
        if (position >= length_)
        {
            throw std::out_of_range("[Get Base of Segment] Position out of range");
        }

        if (!sequence_.empty())
        {
            size_t byte_index = position / (WORD_BIT / 2);
            size_t bit_offset = (position % (WORD_BIT / 2)) * 2;
            uint8_t encoded = (sequence_[byte_index] >> (WORD_BIT - 2 - bit_offset)) & 0x03;
            return GlobalVariant::bin2char[encoded];
        }

        return 'N'; // Default return for empty sequence
    }

    /**
     * @brief Set the Id of the Segment
     * @param {uword_t} id - Segment ID
     */
    void setId(const uword_t id) { id_ = id; }

    /**
     * @brief Set the Length of the Segment
     * @param {uword_t} length - Segment length
     */
    void setSequence(const std::vector<uword_t> &sequence)
    {
        sequence_ = sequence;
        // FIXME: length_ may not be accurate here
        length_ = sequence.size() * (WORD_BIT / 2);
    }

    /**
     * @brief Set the Sequence of the Segment from a string
     * @param {string} &sequence_str - Segment sequence string
     */
    void setSequenceFromString(const std::string &sequence_str)
    {
        length_ = sequence_str.length();
        stringToBinary(sequence_str);
    }

    /**
     * @brief Set the Dis of the Segment
     * @param {uword_t} dis - Segment distance
     */
    void setDis(const uword_t dis) { dis_ = dis; }

    /**
     * @brief Set the SOSR of the Segment
     *  Used only in generating Segment object. It will change the offset and rank info, so be careful when using it.
     * @param {uword_t} offset - Segment offset
     * @param {uword_t} rank - Segment rank
     */
    void setSosr(const uword_t offset, const uword_t rank)
    {
        so_ = offset;
        sr_ = rank;
        has_sosr_ = true;
    }

    void setSnId(const uword_t sn_id)
    {
        sn_id_ = sn_id;
        has_sosr_ = true;
    }

    void setSosrStorage(bool has_sosr, bool compact_mode)
    {
        has_sosr_ = has_sosr;
        sosr_compact_ = compact_mode;
        setStoreSize();
    }

    /**
     * @brief Convert the Segment to GFA text format
     * @return {string} GFA text representation of the segment
     */
    std::string toGfaTxt() const
    {
        std::ostringstream oss;
        oss << "S\t" << id_ << "\t";
        for (size_t i = 0; i < length_; ++i)
        {
            oss << getBaseAt(i);
        }

        oss << "\tSL:i:" << getLength();
        oss << "\tSO:i:" << getOffset();
        oss << "\tSR:i:" << getRank();
        oss << "\tSN:Z:" << getSnId();

        return oss.str();
    }

    /**
     * @brief Convert the Segment to bGFA binary format
     * @param {ofstream} &output_file - Output binary file stream
     * @param {uword_t} &s_offset - Current offset in the binary file
     */
    void toBgfaBin(std::ofstream &output_file, bool store_id = true) const
    {
        if (!output_file.is_open())
        {
            throw std::runtime_error(
                "[Write to bGFA File] Output file is not open");
        }

        uint8_t type_flag = 1; // Type flag for segment

        // output_file.write(reinterpret_cast<const char *>(&type_flag), sizeof(uword_t));
        if (store_id)
        {
            output_file.write(reinterpret_cast<const char *>(&id_),
                              sizeof(uword_t));
        }
        if (has_sosr_)
        {
            if (sosr_compact_)
            {
                if (sr_ > 0x07 || so_ > 0x003FFFFF || sn_id_ > 0x7F)
                {
                    std::cerr << "[BGFA][Warn] SOSR compact overflow: SR(0..7), SO(0..4194303), SN(0..127). "
                              << "Use --sosr-compact=false (uncompressed) to avoid loss." << std::endl;
                }
                const uword_t sr = (sr_ & 0x07);
                const uword_t so = (so_ & 0x003FFFFF);
                const uword_t sn = (sn_id_ & 0x7F);
                const uword_t packed = (sn << 25) | (so << 3) | sr;
                output_file.write(reinterpret_cast<const char *>(&packed),
                                  sizeof(uword_t));
            }
            else
            {
                const uword_t sr = (sr_ & 0x7F);
                const uword_t so = (so_ & 0x01FFFFFF);
                const uword_t packed = (so << 7) | sr;
                output_file.write(reinterpret_cast<const char *>(&packed),
                                  sizeof(uword_t));
                output_file.write(reinterpret_cast<const char *>(&sn_id_),
                                  sizeof(uword_t));
            }
        }

        if (store_id)
        {
            if (length_ > (WORD_BIT - 6) / 2)
            {
                output_file.write(reinterpret_cast<const char *>(&length_),
                                  sizeof(uword_t));

                for (const auto &value : sequence_)
                {
                    output_file.write(reinterpret_cast<const char *>(&value),
                                      sizeof(uword_t));
                }
            }
            else
            {
                uword_t val = (((uword_t(1) << 5) | length_) << (WORD_BIT - 6)) |
                              (sequence_[0] >> 6);
                output_file.write(reinterpret_cast<const char *>(&val),
                                  sizeof(uword_t));
            }
        }
        else
        {
            const uword_t compact1_cap = (WORD_BIT - 6) / 2;
            const uword_t compact1_len_max = std::min<uword_t>(compact1_cap, 0x0F);
            const uword_t compact2_cap = (2 * WORD_BIT - 8) / 2;
            const uword_t compact2_len_max = std::min<uword_t>(compact2_cap, 0x3F);

            if (length_ > compact2_len_max)
            {
                // 00: plain length + sequence words
                output_file.write(reinterpret_cast<const char *>(&length_),
                                  sizeof(uword_t));
                for (const auto &value : sequence_)
                {
                    output_file.write(reinterpret_cast<const char *>(&value),
                                      sizeof(uword_t));
                }
            }
            else if (length_ <= compact1_len_max)
            {
                // 10: [2-bit tag=10][4-bit length][WORD_BIT-6 bits sequence]
                uword_t payload = 0;
                for (uword_t i = 0; i < length_; ++i)
                {
                    const uint8_t enc = GlobalVariant::char2bin[static_cast<uint8_t>(getBaseAt(i))] & 0x3;
                    const uword_t shift = (WORD_BIT - 8) - 2 * i;
                    payload |= (static_cast<uword_t>(enc) << shift);
                }
                const uword_t val = (uword_t(2) << (WORD_BIT - 2)) |
                                    ((length_ & 0x0F) << (WORD_BIT - 6)) |
                                    payload;
                output_file.write(reinterpret_cast<const char *>(&val),
                                  sizeof(uword_t));
            }
            else
            {
                // 11: [2-bit tag=11][6-bit length][WORD_BIT-8 bits sequence] + [next uword sequence]
                uword_t payload_hi = 0;
                uword_t payload_lo = 0;
                const uword_t hi_bits = WORD_BIT - 8;

                for (uword_t i = 0; i < length_; ++i)
                {
                    const uint8_t enc = GlobalVariant::char2bin[static_cast<uint8_t>(getBaseAt(i))] & 0x3;
                    const uword_t bit_pos = 2 * i;

                    const uint8_t b0 = (enc >> 1) & 0x1;
                    const uint8_t b1 = enc & 0x1;

                    const uword_t p0 = bit_pos;
                    const uword_t p1 = bit_pos + 1;

                    if (p0 < hi_bits)
                    {
                        if (b0)
                            payload_hi |= (uword_t(1) << (hi_bits - 1 - p0));
                    }
                    else
                    {
                        const uword_t q = p0 - hi_bits;
                        if (b0)
                            payload_lo |= (uword_t(1) << (WORD_BIT - 1 - q));
                    }

                    if (p1 < hi_bits)
                    {
                        if (b1)
                            payload_hi |= (uword_t(1) << (hi_bits - 1 - p1));
                    }
                    else
                    {
                        const uword_t q = p1 - hi_bits;
                        if (b1)
                            payload_lo |= (uword_t(1) << (WORD_BIT - 1 - q));
                    }
                }

                const uword_t val = (uword_t(3) << (WORD_BIT - 2)) |
                                    ((length_ & 0x3F) << (WORD_BIT - 8)) |
                                    payload_hi;
                output_file.write(reinterpret_cast<const char *>(&val),
                                  sizeof(uword_t));
                output_file.write(reinterpret_cast<const char *>(&payload_lo),
                                  sizeof(uword_t));
            }
        }
    }

    /**
     * @brief Equality operator for Segment
     * @details Two segments are equal if their sequences are same.
     * @param {Segment} &other - Other segment to compare
     * @return {bool} True if segments are equal, false otherwise
     */
    bool operator==(const Segment &other) const
    {
        return sequence_ == other.sequence_;
    }

    /**
     * @brief Inequality operator for Segment
     * @param {Segment} &other - Other segment to compare
     * @return {bool} True if segments are not equal, false otherwise
     */
    bool operator!=(const Segment &other) const
    {
        return !(*this == other);
    }

    /**
     * @brief Less-than operator for Segment (based on distance, then ID)
     * @param {Segment} &other - Other segment to compare
     * @return {bool}
     */
    bool operator<(const Segment &other) const
    {
        if (dis_ != other.dis_)
        {
            return dis_ < other.dis_;
        }
        return id_ < other.id_;
    }

    /**
     * @brief Greater-than operator for Segment
     * @param {Segment} &other - Other segment to compare
     * @return {bool}
     */
    bool operator>(const Segment &other) const { return other < *this; }

    /**
     * @brief Less-than-or-equal-to operator for Segment
     * @param {Segment} &other - Other segment to compare
     * @return {bool}
     */
    bool operator<=(const Segment &other) const { return !(other < *this); }

    /**
     * @brief Greater-than-or-equal-to operator for Segment
     * @param {Segment} &other - Other segment to compare
     * @return {bool}
     */
    bool operator>=(const Segment &other) const { return !(*this < other); }

    /**
     * @brief Output stream operator for Segment
     * @param {ostream} &os - Output stream
     * @param {Segment} &segment - Segment to output
     * @return {ostream&} Reference to the output stream
     */
    friend std::ostream &operator<<(std::ostream &os, const Segment &segment)
    {
        // std::ostringstream oss;
        os << "S: " << segment.getId()
           << ", Lengsth: " << segment.getLength();

        size_t preview_len = std::min(uword_t(100), segment.getLength());
        os << ", Preview: ";
        for (size_t i = 0; i < preview_len; ++i)
        {
            os << segment.getBaseAt(i);
        }
        if (segment.getLength() > preview_len)
        {
            os << "...";
        }

        os << "\n\tSO: " << segment.getOffset()
           << "\tSR: " << segment.getRank();

        if (segment.dis_ > 0)
        {
            os << "\tDis: " << segment.dis_;
        }

        os << std::endl;

        return os;
    }

private:
    /**
     * @brief Set the Sequence of the Segment from a string
     * @param {string} &seq_str - Segment sequence string
     */
    void stringToBinary(const std::string &seq_str)
    {
        uint8_t max_number = WORD_BIT / 2;
        size_t num_bytes = (seq_str.length() + max_number - 1) / max_number;

        sequence_.resize(num_bytes, 0);
        for (size_t i = 0; i < seq_str.length(); ++i)
        {
            uint8_t encoded = GlobalVariant::char2bin[seq_str[i]];
            size_t byte_index = i / max_number;
            size_t bit_offset = (i % max_number) * 2;
            sequence_[byte_index] |= (encoded << (WORD_BIT - 2 - bit_offset));
        }
    }

    void setStoreSize()
    {
        store_size_ = 1; // for id

        if (has_sosr_)
            store_size_ += sosr_compact_ ? 1 : 2; // sosr (+sn)

        if (length_ > (WORD_BIT - 6) / 2)
        {
            store_size_ += 1; // for length
            store_size_ += sequence_.size();
        }
        else
        {
            store_size_ += 1; // compact representation
        }
    }
};

/**
 * @brief LinksSet class representing a set of links in the graph
 */
class LinksSet
{
private:
    uword_t link_num_;                               ///< Number of links
    uword_t store_size_;                             ///< Size in storage (in uword_t units)
    std::vector<std::unordered_set<uword_t>> links_; ///< Adjacency set representation of links, links[from_id] contains encoded (to_id<<2)|dir

public:
    /**
     * @brief Default constructor initializes an empty LinksSet
     */
    LinksSet() : link_num_(0), store_size_(0) {}

    /**
     * @brief Constructor that builds LinksSet from given data and indptr vectors
     * @param {uword_t} links_num - Number of links
     * @param {uword_t} nodes_num - Number of nodes
     * @param {vector<uword_t>} &data - Data vector containing link information
     * @param {vector<uword_t>} &indptr - Indptr vector for indexing
     */
    explicit LinksSet(const uword_t links_num, const uword_t nodes_num,
                      const std::vector<uword_t> &data,
                      const std::vector<uword_t> &indptr)
        : link_num_(links_num)
    {
        links_.resize(nodes_num);
        for (uword_t i = 0; i < nodes_num; ++i)
        {
            uword_t start = indptr[i];
            uword_t end = indptr[i + 1];
            for (uword_t j = start; j < end; ++j)
            {
                links_[i].insert(data[j]);
            }
        }
        setStoreSize();
    }

    /**
     * @brief Get the number of links
     * @return {uword_t} Number of links
     */
    uword_t getLinkNum() const { return link_num_; }

    /**
     * @brief Get the links from a specific node
     * @param {uword_t} index - Node index
     * @return {const vector<uword_t>&} Vector of links from the specified node
     */
    const std::unordered_set<uword_t> &getLinksFromNode(uword_t index) const
    {
        if (index >= links_.size())
        {
            throw std::out_of_range("[Get Link Info] Node index out of range");
        }
        return links_[index];
    }

    /**
     * @brief Get all links in the LinksSet
     * @return {const vector<vector<uword_t>>&} 2D vector of all links
     */
    const std::vector<std::unordered_set<uword_t>> &getAllLinks() const
    {
        return links_;
    }

    /**
     * @brief Get the target node index of a link by link index
     * @param {uword_t} from_id - Source node index
     * @param {uword_t} i - Link index
     * @return {uword_t} Target node index
     */
    uword_t getLinkToByIndex(uword_t from_id, uword_t i) const
    {
        if (from_id >= links_.size())
        {
            throw std::out_of_range("[Get Link Info] Node index out of range");
        }
        if (i >= links_[from_id].size())
        {
            throw std::out_of_range("[Get Link Info] Link index out of range");
        }
        auto it = links_[from_id].begin();
        std::advance(it, i);
        return ((*it) >> 2);
    }

    /**
     * @brief Get the direction of a link from a specific node to a target node
     * @param {uword_t} from_id - Source node index
     * @param {uword_t} to_id - Target node index
     * @return {uword_t} Direction of the link (2-bit encoded)
     */
    uword_t getLinkDir(uword_t from_id, uword_t to_id) const
    {
        if (from_id >= links_.size())
        {
            throw std::out_of_range("[Get Link Info] Node index out of range");
        }
        for (const auto &link : links_[from_id])
        {
            if ((link >> LINK_DIR_BIT) == to_id)
            {
                return (link & 0x03);
            }
        }
        throw std::logic_error("[Get Link Info] Link not found");
    }

    /**
     * @brief Get the direction of a link from a specific node by link index
     * @param {uword_t} from_id - Source node index
     * @param {uword_t} i - Link index
     * @return {uword_t} Direction of the link (2-bit encoded)
     */
    uword_t getLinkDirByIndex(uword_t from_id, uword_t i) const
    {
        if (from_id >= links_.size())
        {
            throw std::out_of_range("[Get Link Info] Node index out of range");
        }
        if (i >= links_[from_id].size())
        {
            throw std::out_of_range("[Get Link Info] Link index out of range");
        }
        auto it = links_[from_id].begin();
        std::advance(it, i);
        return ((*it) & 0x03);
    }

    /**
     * @brief Get the storage size of the LinksSet in uword_t units
     * @return {uword_t} Storage size
     */
    uword_t getStoreSize()
    {
        setStoreSize();
        return store_size_;
    }

    /**
     * @brief Add a new empty node to the LinksSet
     */
    void addNode()
    {
        links_.emplace_back();
        setStoreSize();
    }

    /**
     * @brief Resize the LinksSet to accommodate a new size
     * @param {uword_t} new_size - New size for the LinksSet
     */
    void resizeLinks(uword_t new_size)
    {
        if (new_size >= links_.size())
        {
            links_.resize(new_size);
            setStoreSize();
        }
    }

    /**
     * @brief Add a link from from_id to to_id with specified directions
     * @param {uword_t} from_id - Source node index
     * @param {uword_t} to_id - Target node index
     * @param {string} from_dir - Direction from source node ('+' or '-')
     * @param {string} to_dir - Direction to target node ('+' or '-')
     */
    void addLink(uword_t from_id, uword_t to_id,
                 std::string from_dir, std::string to_dir)
    {
        if (from_id >= links_.size() || to_id >= links_.size())
        {
            throw std::out_of_range("[Add Link] Node index out of range");
        }
        if (from_dir.length() != 1 || to_dir.length() != 1 ||
            (from_dir[0] != '+' && from_dir[0] != '-') ||
            (to_dir[0] != '+' && to_dir[0] != '-'))
        {
            throw std::invalid_argument("[Add Link] Invalid direction format");
        }
        uword_t to_value = (GlobalVariant::dir2bin[from_dir[0]] << 1) | GlobalVariant::dir2bin[to_dir[0]];
        to_value |= (to_id << 2);
        auto res = links_[from_id].insert(to_value);
        if (res.second)
            link_num_++;

        setStoreSize();
    }

    /**
     * @brief Add a link from from_id to to_value
     * @param {uword_t} from_id - Source node index
     * @param {uword_t} to_value - Encoded target node index and direction
     */
    void addLink(uword_t from_id, uword_t to_value)
    {
        if (from_id >= links_.size() || (to_value >> 2) >= links_.size())
        {
            throw std::out_of_range("[Add Link] Node index out of range");
        }
        auto res = links_[from_id].insert(to_value);
        if (res.second)
            link_num_++;

        setStoreSize();
    }

    /**
     * @brief Redirect all incoming links pointing to old_to so they point to new_to.
     * @param old_to - Original target node id
     * @param new_to - New target node id
     */
    void redirectIncoming(uword_t old_to, uword_t new_to)
    {
        if (old_to >= links_.size() || new_to >= links_.size())
        {
            throw std::out_of_range("[Redirect Incoming] Node index out of range");
        }
        for (auto &from_links : links_)
        {
            std::vector<uword_t> to_replace;
            for (const auto &to_value : from_links)
            {
                uword_t to = (to_value >> LINK_DIR_BIT);
                if (to == old_to)
                {
                    to_replace.push_back(to_value);
                }
            }
            for (const auto &old_val : to_replace)
            {
                uword_t dir_bits = (old_val & 0x03);
                uword_t new_val = (new_to << LINK_DIR_BIT) | dir_bits;
                // erase old
                size_t erased = from_links.erase(old_val);
                if (erased > 0)
                {
                    link_num_--; // one edge removed
                }
                // insert new (may or may not exist)
                auto res = from_links.insert(new_val);
                if (res.second)
                {
                    link_num_++; // new unique edge added
                }
            }
        }

        setStoreSize();
    }

    /**
     * @brief Move all outgoing links from old_from to new_from, preserving direction bits.
     *        After moving, the old_from's outgoing links are cleared. The total link count remains unchanged.
     * @param old_from - Original source node id
     * @param new_from - New source node id
     */
    void redirectOutgoing(uword_t old_from, uword_t new_from)
    {
        if (old_from >= links_.size() || new_from >= links_.size())
        {
            throw std::out_of_range("[Redirect Outgoing] Node index out of range");
        }
        if (old_from == new_from)
            return;
        auto &src = links_[old_from];
        auto &dst = links_[new_from];
        uword_t inserted = 0;
        for (const auto &val : src)
        {
            auto res = dst.insert(val);
            if (res.second)
                inserted++;
        }
        // remove all from src
        uword_t removed = src.size();
        src.clear();
        // update unique link count: +inserted (new unique) - removed (old unique)
        if (removed > inserted)
            link_num_ -= (removed - inserted);
        else if (inserted > removed)
            link_num_ += (inserted - removed);

        setStoreSize();
    }

    /**
     * @brief Remove all outgoing links from a node.
     * @param from_id - Source node id
     */
    void removeOutgoing(uword_t from_id)
    {
        if (from_id >= links_.size())
        {
            throw std::out_of_range("[Remove Outgoing] Node index out of range");
        }
        link_num_ -= links_[from_id].size();
        links_[from_id].clear();

        setStoreSize();
    }

    /**
     * @brief Remove all incoming links to a node.
     * @param to_id - Target node id
     */
    void removeIncoming(uword_t to_id)
    {
        if (to_id >= links_.size())
        {
            throw std::out_of_range("[Remove Incoming] Node index out of range");
        }
        for (auto &from_links : links_)
        {
            for (auto it = from_links.begin(); it != from_links.end();)
            {
                if (((*it) >> LINK_DIR_BIT) == to_id)
                {
                    it = from_links.erase(it);
                    link_num_--;
                }
                else
                {
                    ++it;
                }
            }
        }

        setStoreSize();
    }

    /**
     * @brief Convert the LinksSet to GFA text format
     * @return {string} GFA text representation of the links
     */
    std::string toGfaTxt() const
    {
        std::ostringstream oss;
        uword_t link_id = 0;
        for (size_t i = 0; i < links_.size(); ++i)
        {
            for (const auto &link : links_[i])
            {
                oss << "L\t" << i << "\t"
                    << GlobalVariant::bin2dir[(link & 0x03) >> 1] << "\t"
                    << (link >> 2) << "\t"
                    << GlobalVariant::bin2dir[link & 0x01]
                    << "\t0M" << std::endl;
                // link_id++;
            }
        }
        return oss.str();
    }

    /**
     * @brief Convert the LinksSet to bGFA binary format
     * @param {ofstream} &output_file - Output binary file stream
     * @param {uword_t} &l_offset - Current offset in the binary file
     */
    void toBgfaBin(std::ofstream &output_file) const
    {
        if (!output_file.is_open())
        {
            throw std::runtime_error("[Write to bGFA File] Output file is not open");
        }

        uword_t n_rows = links_.size();

#ifdef INFO_DEBUG
        std::cout << "[LinksSet to bGFA Bin] link_num: " << link_num_
                  << ", n_rows(links_.size()): " << n_rows
                  << ", store_size: " << store_size_ << std::endl;
#endif
        // size of links
        output_file.write(reinterpret_cast<const char *>(&link_num_),
                          sizeof(uword_t));
        // size of nodes(index)
        output_file.write(reinterpret_cast<const char *>(&n_rows),
                          sizeof(uword_t));

        uword_t indptr = 0;
        output_file.write(reinterpret_cast<const char *>(&indptr),
                          sizeof(uword_t));

        for (size_t i = 0; i < n_rows; ++i)
        {
            indptr = indptr + links_[i].size();
            output_file.write(reinterpret_cast<const char *>(&indptr),
                              sizeof(uword_t));
        }

        for (size_t i = 0; i < n_rows; ++i)
        {
            for (const auto &link : links_[i])
            {
                output_file.write(reinterpret_cast<const char *>(&link),
                                  sizeof(uword_t));
                // current_ptr++;
            }
            // indptr[i + 1] = current_ptr;
        }
    }

    /**
     * @brief Output stream operator for LinksSet
     * @param {ostream} &os - Output stream
     * @param {LinksSet} &links_set - LinksSet to output
     * @return {ostream&} Reference to the output stream
     */
    friend std::ostream &operator<<(std::ostream &os, const LinksSet &links_set)
    {
        uword_t link_id = 0;
        for (size_t i = 0; i < links_set.links_.size(); ++i)
        {
            for (const auto &link : links_set.links_[i])
            {
                os << "L " << i << " "
                   << GlobalVariant::bin2dir[(link & 0x03) >> 1] << " "
                   << (link >> 2) << " "
                   << GlobalVariant::bin2dir[link & 0x01] << std::endl;
                // link_id++;
            }
        }
        return os;
    }

private:
    void setStoreSize()
    {
        store_size_ = 2;                  // for link_num_ and n_rows
        store_size_ += links_.size() + 1; // for indptr
        store_size_ += link_num_;         // for data
    }
};

class Path
{
private:
    uword_t id_;
    /**
     * @brief Nodes in the Path
     * @details If Direct Mode, stored as encoded (node_id << 1) | dir
     *          If Increment Mode: flag_seg, head_seg, inc_seg, inc_seg, ...
     *              flag: [31 bits unused(max 15 dirs)][1bit flag]
     *                  2bits->dir: 00:终止, 10:-, 11:+
     *                  1bit->flag: 0: seg_mode, 1: inc_mode
     *              head_seg: encoded (node_id << 1) | dir
     *              inc_seg:  4 x inc(8bits, max 255)
     */
    std::vector<uword_t> nodes_;
    uword_t mode_;
    uword_t store_size_; ///< Size in storage (in uword_t units)

public:
    Path() = default;
    explicit Path(const std::vector<uword_t> &nodes) : nodes_(nodes), mode_(PATHMODE_DIRECT)
    {
        setStoreSize();
    }

    explicit Path(uword_t id, const std::vector<uword_t> &nodes) : id_(id), nodes_(nodes), mode_(PATHMODE_DIRECT)
    {
        setStoreSize();
    }

    explicit Path(uword_t path_id, std::string path_str, NodeMap &id_to_id, PathMode mode) : id_(path_id), mode_(mode)
    {
        if (mode_ == PATHMODE_DIRECT)
            parsePathDirect(path_str, id_to_id);
        else if (mode_ == PATHMODE_INCREMENT)
        {
            std::vector<std::pair<uword_t, char>> tokens;
            {
                std::istringstream iss(path_str);
                std::string node_str;
                while (std::getline(iss, node_str, ','))
                {
                    if (node_str.empty())
                    {
                        throw std::logic_error(
                            "[Create Paths Info] Invalid path string: " + path_str);
                    }
                    char dir = node_str.back();
                    uword_t node_id = id_to_id[node_str.substr(0, node_str.length() - 1)];
                    tokens.emplace_back(node_id, dir);
                }
            }

            parsePathIncrement(tokens);
        }
    }

    explicit Path(std::string walk_str, NodeMap &id_to_id, PathMode mode) : id_(0), mode_(mode)
    {
        // Parse walk_str into tokens: (>|<)id(>|<)id...(>|<)id
        auto is_space = [](char c) -> bool
        { return c == ' ' || c == '\t' || c == '\r' || c == '\n'; };

        size_t pos = 0;
        std::vector<std::pair<uword_t, char>> tokens;
        auto read_next = [&]() -> bool
        {
            // skip whitespaces
            while (pos < walk_str.size() && is_space(walk_str[pos]))
                ++pos;
            if (pos >= walk_str.size())
                return false;

            char dir_sym = walk_str[pos];
            if (dir_sym != '>' && dir_sym != '<')
            {
                throw std::logic_error(
                    "[Create Paths Info] Invalid walk token direction in: " + walk_str);
            }
            ++pos;

            // read id until next direction symbol or whitespace
            size_t start = pos;
            while (pos < walk_str.size() && !is_space(walk_str[pos]) && walk_str[pos] != '>' && walk_str[pos] != '<')
                ++pos;
            if (start == pos)
            {
                throw std::logic_error(
                    "[Create Paths Info] Missing id after direction in walk: " + walk_str);
            }
            std::string id_str = walk_str.substr(start, pos - start);

            uword_t node_id = id_to_id[id_str];
            char dir_char = (dir_sym == '>') ? '+' : '-';
            tokens.emplace_back(node_id, dir_char);
            return true;
        };

        while (read_next())
            ;

        if (tokens.empty())
        {
            setStoreSize();
            return;
        }

        if (mode_ == PATHMODE_DIRECT)
        {
            // Build direct-mode nodes: encoded (id<<1)|dir_bit
            uword_t first = packIdWithDir(tokens[0].first, tokens[0].second);
            nodes_.push_back(first); // from
            nodes_.push_back(first); // placeholder end
            for (size_t i = 1; i < tokens.size(); ++i)
            {
                uword_t node = packIdWithDir(tokens[i].first, tokens[i].second);
                nodes_.push_back(node);
            }
            nodes_[1] = nodes_.back();
        }
        else if (mode_ == PATHMODE_INCREMENT)
        {
            // Increment mode: group into blocks and pack dirs/increments
            parsePathIncrement(tokens);
        }
        else
        {
            throw std::logic_error("[Create Paths Info] Unknown path mode for walk parsing");
        }

        setStoreSize();
    }

    // uword_t getLength() const { return nodes_.size(); }

    /**
     * @brief Get the from node of the Path
     * @return {uword_t} From node
     */
    uword_t getFromNode() const
    {
        return nodes_[0];
    }

    /**
     * @brief Get the end node of the Path
     * @return {uword_t} End node
     */
    uword_t getEndNode() const { return nodes_[1]; }

    /**
     * @brief Get all nodes in the Path
     * @return {const vector<uword_t>&} Vector of all nodes
     */
    const std::vector<uword_t> &getAllNodes() const { return nodes_; }

    uword_t getId() const { return id_; }

    PathMode getMode() const { return static_cast<PathMode>(mode_); }

    /**
     * @brief Get the storage size of the Path in uword_t units
     * @return {uword_t} Storage size
     */
    uword_t getStoreSize()
    {
        setStoreSize();
        return store_size_;
    }

    void setMode(PathMode mode)
    {
        if (mode != PATHMODE_DIRECT && mode != PATHMODE_INCREMENT)
        {
            throw std::invalid_argument("[Set Path Mode] Invalid path mode");
        }
        mode_ = mode;
        setStoreSize();
    }

    /**
     * @brief Convert the Path to GFA text format
     * @param {LinksSet} &links_set - LinksSet for reference
     * @return {string} GFA text representation of the path
     */
    std::string toGfaTxt() const
    {
        if (mode_ == PATHMODE_DIRECT)
        {
            return toGfaTxtDirect();
        }
        else if (mode_ == PATHMODE_INCREMENT)
        {
            return toGfaTxtIncrement();
        }
        else
        {
            throw std::logic_error("[Path to GFA Txt] Unknown path mode");
        }
    }

    /**
     * @brief Convert the Path to bGFA binary format
     * @param {ofstream} &output_file - Output binary file stream
     */
    void toBgfaBin(std::ofstream &output_file) const
    {
        if (!output_file.is_open())
        {
            throw std::runtime_error(
                "[Write to bGFA File] Output file is not open");
        }

        uword_t n_nodes = nodes_.size();

        output_file.write(reinterpret_cast<const char *>(&n_nodes),
                          sizeof(uword_t));

        for (const auto &node : nodes_)
        {
            output_file.write(reinterpret_cast<const char *>(&node),
                              sizeof(uword_t));
        }
    }

    /**
     * @brief Output stream operator for Path
     * @param {ostream} &os - Output stream
     * @param {Path} &path - Path to output
     * @return {ostream&} Reference to the output stream
     */
    friend std::ostream &operator<<(std::ostream &os, const Path &path)
    {
        os << "P";

        if (path.mode_ == PATHMODE_DIRECT)
        {
            os << "\t(Direct Mode)\t";
            os << ", From Node: " << (path.nodes_[0] >> 1)
               << ", End Node: " << (path.nodes_[1] >> 1)
               << "\n\t" << (path.nodes_[0] >> 1)
               << GlobalVariant::bin2dir[path.nodes_[0] & 0x01] << " -> ";

            for (size_t i = 2; i < path.nodes_.size(); ++i)
            {
                os << (path.nodes_[i] >> 1)
                   << GlobalVariant::bin2dir[path.nodes_[i] & 0x01];
                if (i != path.nodes_.size() - 1)
                    os << " -> ";
            }
        }
        else if (path.mode_ == PATHMODE_INCREMENT)
        {
            os << "\t(Increment Mode)\t";
            os << ", From Node: " << (path.nodes_[0] >> 2)
               << ", End Node: " << (path.nodes_[1] >> 2) << "\n";

            os << path.toGfaTxtIncrement();
        }
        return os;
    }

private:
    /**
     * @brief Check if the link from from_id to to_id with specified directions is not unique
     * @param {LinksSet} &links_set - LinksSet to check
     * @param {uword_t} from_id - Source node index
     * @param {uword_t} to_id - Target node index
     * @param {char} from_dir - Direction from source node ('+' or '-')
     * @param {char} to_dir - Direction to target node ('+' or '-')
     * @return {bool} True if the link is not unique, false otherwise
     */
    bool isNotUniqueSuccessor(LinksSet &links_set, uword_t from_id, uword_t to_id, char from_dir, char to_dir) const
    {
        uword_t expected_dir = (GlobalVariant::dir2bin[from_dir] << 1) |
                               (GlobalVariant::dir2bin[to_dir]);

        if (from_id >= links_set.getAllLinks().size())
        {
            throw std::out_of_range(
                "[Check Is Unique Successor] Node index out of range");
        }

        const auto &links = links_set.getLinksFromNode(from_id);
        if (links.size() == 1)
        {
            if (links_set.getLinkToByIndex(from_id, 0) == to_id &&
                links_set.getLinkDirByIndex(from_id, 0) == expected_dir)
                return false;
            else
                throw std::logic_error(
                    "[Check Is Unique Successor] Link not found" +
                    std::to_string(from_id) + from_dir + " -> " +
                    std::to_string(to_id) + to_dir);
        }
        for (int i = 0; i < links.size(); ++i)
        {
            if (links_set.getLinkToByIndex(from_id, i) == to_id &&
                links_set.getLinkDirByIndex(from_id, i) == expected_dir)
                return true;
        }

        throw std::logic_error(
            "[Check Is Unique Successor] Link not found: " +
            std::to_string(from_id) + from_dir + " -> " +
            std::to_string(to_id) + to_dir);
    }

    void parsePathDirect(std::string path_str, NodeMap &id_to_id)
    {
        std::istringstream iss(path_str);
        std::string node_str;
        uword_t node, node_id;
        char dir;

        // from id
        std::getline(iss, node_str, ',');
        if (node_str.empty())
        {
            throw std::logic_error(
                "[Create Paths Info] Invalid path string: " + path_str);
        }
        dir = node_str.back();
        node_id = id_to_id[node_str.substr(0, node_str.length() - 1)];
        node = packIdWithDir(node_id, dir);
        nodes_.push_back(node);
        nodes_.push_back(node);

        // to ids
        std::getline(iss, node_str, ',');

        while (1)
        {
            // next node
            dir = node_str.back();
            node_id = id_to_id[node_str.substr(0, node_str.length() - 1)];
            node = packIdWithDir(node_id, dir);

            nodes_.push_back(node);

            if (!std::getline(iss, node_str, ','))
            {
                nodes_[1] = node;
                break;
            }
        }
    }

    void parsePathIncrement(std::vector<std::pair<uword_t, char>> tokens)
    {
        uword_t idx = 0;
        uword_t node;

        node = packIdDirFlag(tokens[idx].first, tokens[idx].second, 0);
        nodes_.push_back(node); // nodes[0], seg[0]
        node = packIdDirFlag(tokens[tokens.size() - 1].first, tokens[tokens.size() - 1].second, 0);
        nodes_.push_back(node); // nodes[1], seg[-1] (end placeholder)
        idx++;                  // idx = 1

        // If only one node, finalize
        if (tokens.size() == 1)
        {
            return;
        }

        // Save second node at nodes_[2]
        node = packIdDirFlag(tokens[idx].first, tokens[idx].second, 0);
        nodes_.push_back(node); // nodes[2], seg1
        idx++;                  // idx = 2

        uword_t inc_head_pos = idx;           // position of "head" before a potential increment group
        uword_t inc_pre_id = tokens[1].first; // previous id to compute deltas

        // State for an ongoing increment block
        bool in_block = false;
        uword_t inc_word = 0;       // packs 4 x 8-bit increments
        uword_t inc_count = 0;      // number of increments recorded in current block (max 16)
        uword_t dir_word_index = 0; // nodes_ index where direction word is stored when in_block

        while (idx + 1 < tokens.size())
        {
            // Consider starting from tokens[idx]
            uword_t cur_id = tokens[idx].first;
            char cur_dir = tokens[idx].second;
            uword_t delta = cur_id < inc_pre_id ? PATH_MAX_INC + 1 : (cur_id - inc_pre_id);

            if (!in_block)
            {
                // 1.完整模式, 保存当前段，continue
                if (delta > PATH_MAX_INC)
                {
                    uword_t node = packIdDirFlag(cur_id, cur_dir, 0);
                    nodes_.push_back(node);
                    inc_head_pos = nodes_.size() - 1;
                    inc_pre_id = cur_id;
                    idx++;
                    continue;
                }

                // 2.可以采用增量，但是不具有4个可用段（含当前）, 保存所有段，break
                if (idx + PATH_INC_SIZE - 1 >= tokens.size())
                {
                    for (; idx < tokens.size(); ++idx)
                    {
                        uword_t node = packIdDirFlag(tokens[idx].first, tokens[idx].second, 0);
                        nodes_.push_back(node);
                        inc_head_pos = nodes_.size() - 1;
                        inc_pre_id = tokens[idx].first;
                    }
                    break;
                }

                // 3.可以采用增量，具有4个可用段（含当前）
                uword_t prev = inc_pre_id;
                bool ok = true;
                inc_word = 0;
                uword_t dir_word = 0;
                for (int k = 0; k < PATH_INC_SIZE; ++k)
                {
                    uword_t nid = tokens[idx + k].first;
                    if (nid > PATH_MAX_INC + prev)
                    {
                        ok = false;
                        inc_word = 0;
                        break;
                    }
                    dir_word = setDirWord(dir_word, tokens[idx + k].second, k);
                    inc_word = setIncWord(inc_word, nid - prev, k);
                    prev = nid;
                }

                // 3.1.四个段增量超出范围，保存当前段，continue
                if (!ok)
                {
                    uword_t node = packIdDirFlag(cur_id, cur_dir, 0);
                    nodes_.push_back(node);
                    inc_head_pos = nodes_.size() - 1;
                    inc_pre_id = cur_id;
                    idx++;
                    continue;
                }

                // 3.2.四个段增量未超出范围
                // 为增量开始头结点进行标记
                nodes_[inc_head_pos] |= 1u;

                // 打包四个dir到dir_word
                nodes_.push_back(dir_word);
                dir_word_index = nodes_.size() - 1;

                // 打包四个inc到inc_word
                nodes_.push_back(inc_word);
                inc_word = 0;

                inc_pre_id = prev;
                inc_count = PATH_INC_SIZE;
                in_block = true; // reset for subsequent chunks
                idx += PATH_INC_SIZE;
                continue;
            }
            else
            {
                // Extend current block while delta within PATH_MAX_INC and count <= 16
                if (delta <= PATH_MAX_INC && inc_count < 16)
                {
                    // Append direction (2 bits) into existing direction word
                    nodes_[dir_word_index] =
                        setDirWord(nodes_[dir_word_index], cur_dir, inc_count);

                    // Append increment into the current 4-slot word
                    uword_t slot = inc_count % 4;
                    inc_word = setIncWord(inc_word, delta, slot);
                    inc_count++;
                    inc_pre_id = cur_id;

                    // Every 4 increments, flush inc_word
                    if ((inc_count % 4) == 0)
                    {
                        nodes_.push_back(inc_word);
                        inc_word = 0;
                    }
                    idx++;
                    continue;
                }

                // 不满足距离需求，需要将现有块存储（若inc_count为4的倍数已在之前flush，无需再存储）
                if (inc_count % 4)
                {
                    nodes_.push_back(inc_word);
                    inc_word = 0;
                }

                inc_count = 0;
                in_block = false;

                // 存储下一个节点为完整模式
                uword_t node = packIdDirFlag(cur_id, cur_dir, 0);
                nodes_.push_back(node);
                inc_head_pos = nodes_.size() - 1;
                inc_pre_id = cur_id;
                idx++;
                continue;
            }
        }

        // Finalize: flush any pending increments if in a block
        if (in_block && (inc_count % 4))
        {
            nodes_.push_back(inc_word);
        }
    }

    std::string toGfaTxtDirect() const
    {
        std::ostringstream oss;

        oss << "P\t" << id_ << "\t";

        // Print from node
        oss << (nodes_[0] >> 1) << GlobalVariant::bin2dir[nodes_[0] & 0x01];

        // Print intermediate nodes (skip nodes_[1] which is end placeholder)
        for (size_t i = 2; i < nodes_.size(); ++i)
        {
            oss << "," << (nodes_[i] >> 1) << GlobalVariant::bin2dir[nodes_[i] & 0x01];
        }

        return oss.str();
    }

    std::string toGfaTxtIncrement() const
    {
        std::ostringstream oss;

        oss << "P\t" << id_ << "\t";

        if (nodes_.empty())
            return std::string();

        // Print start node (nodes_[0])
        uword_t start_id = (nodes_[0] >> 2);
        char start_dir = GlobalVariant::bin2dir[(nodes_[0] >> 1) & 0x01];
        oss << start_id << start_dir;
        uword_t prev_id = start_id;

        // Iterate over encoded sequence, skipping nodes_[1] (end placeholder)
        size_t i = 2;
        while (i < nodes_.size())
        {
            uword_t val = nodes_[i];
            if ((val & 0x01u) == 1u)
            {
                // Head of an increment block
                uword_t head_id = (val >> 2);
                char head_dir = GlobalVariant::bin2dir[(val >> 1) & 0x01];
                oss << "," << head_id << head_dir;
                prev_id = head_id;

                if (i + 1 >= nodes_.size())
                    break; // malformed; no direction word
                uword_t dir_word = nodes_[i + 1];

                // Count directions encoded: 2 bits per node; 00 indicates termination
                int dir_count = 0;
                for (int p = 0; p < 16; ++p)
                {
                    uword_t code = (dir_word >> (2 * p)) & 0x03u;
                    if (code == 0u)
                        break;
                    ++dir_count;
                }
                int inc_words = (dir_count + PATH_INC_SIZE - 1) / PATH_INC_SIZE;

                // Collect increment words following direction word
                std::vector<uword_t> incs;
                incs.reserve(inc_words);
                size_t base = i + 2;
                for (int w = 0; w < inc_words && (base + w) < nodes_.size(); ++w)
                {
                    incs.push_back(nodes_[base + w]);
                }

                // Reconstruct nodes from increments and directions
                for (int m = 0; m < dir_count; ++m)
                {
                    uword_t inc_word = incs[m / PATH_INC_SIZE];
                    uword_t delta = (inc_word >> (PATH_INC_SHIFT * (m % PATH_INC_SIZE))) & 0xFFu;
                    prev_id += delta;

                    uword_t code = (dir_word >> (2 * m)) & 0x03u; // 2:'-', 3:'+'
                    char dir_char = (code == 3u) ? '+' : '-';
                    oss << "," << prev_id << dir_char;
                }

                // Advance past head, dir_word and inc_words
                i = base + inc_words;
                continue;
            }
            else
            {
                // Normal node stored as full id+dir
                uword_t id = (val >> 2);
                char dir_c = GlobalVariant::bin2dir[(val >> 1) & 0x01];
                oss << "," << id << dir_c;
                prev_id = id;
                ++i;
            }
        }

        return oss.str();
    }

    uword_t packIdWithDir(uword_t id, char dir)
    {
        uword_t res = (id << 1) | (GlobalVariant::dir2bin[dir]);
        return res;
    }

    uword_t packIdDirFlag(uword_t id, char dir, uword_t flag)
    {
        uword_t res = (id << 2) | (GlobalVariant::dir2bin[dir] << 1) | flag;
        return res;
    }

    uword_t setDirWord(uword_t word, char dir, uword_t pos)
    {
        uword_t code = (dir == '+') ? 3u : 2u;
        word |= (code << (2 * pos));
        return word;
    }

    uword_t setIncWord(uword_t word, uword_t inc, uword_t slot)
    {
        word |= ((inc & 0xFFu) << (PATH_INC_SHIFT * slot));
        return word;
    }

    /**
     * @brief Set the storage size of the Path in uword_t units
     */
    void setStoreSize()
    {
        store_size_ = 1;              // for n_nodes
        store_size_ += nodes_.size(); // for nodes
    }
};

class Walk
{
private:
    Path path_;
    uword_t sample_id_;    ///< Sample ID
    uword_t haplotype_id_; ///< Haplotype ID
    uword_t sequence_id_;  ///< Sequence ID
    uword_t store_size_;   ///< Size in storage (in uword_t units)

public:
    Walk() : sample_id_(0), haplotype_id_(0), sequence_id_(0), store_size_(0)
    {
        setStoreSize();
    }

    explicit Walk(const Path &path)
        : path_(path), sample_id_(0), haplotype_id_(0), sequence_id_(0), store_size_(0)
    {
        setStoreSize();
    }

    explicit Walk(uword_t sample_id, uword_t haplotype_id, uword_t sequence_id, const Path &path)
        : path_(path), sample_id_(sample_id), haplotype_id_(haplotype_id), sequence_id_(sequence_id), store_size_(0)
    {
        setStoreSize();
    }

    /**
     * @brief Get the sample ID of the Walk
     * @return {uword_t} Sample ID
     */
    uword_t getSampleId() const { return sample_id_; }

    /**
     * @brief Get the haplotype ID of the Walk
     * @return {uword_t} Haplotype ID
     */
    uword_t getHaplotypeId() const { return haplotype_id_; }

    uword_t getSequenceId() const { return sequence_id_; }

    /**
     * @brief Get the underlying Path of the Walk
     * @return {const Path&} Path object
     */
    const Path &getPath() const { return path_; }

    /**
     * @brief Get the storage size of the Walk in uword_t units
     * @return {uword_t} Storage size
     */
    uword_t getStoreSize()
    {
        setStoreSize();
        return store_size_;
    }

    /**
     * @brief Convert the Walk to GFA text format
     * @param {LinksSet} &links_set - LinksSet for reference (unused)
     * @return {string} GFA text representation of the walk
     */
    std::string toGfaTxt() const
    {
        std::ostringstream oss;
        // Header: sample and haplotype
        oss << "W\t" << sample_id_ << "\t" << haplotype_id_ << "\t";

        // Use Path's textual nodes representation (handles direct/increment)
        std::string nodes_txt = path_.toGfaTxt();
        if (nodes_txt.empty())
            return oss.str();

        // tokens like "id+" or "id-" separated by commas
        size_t pos = 0;
        while (pos < nodes_txt.size())
        {
            size_t next = nodes_txt.find(',', pos);
            std::string tok = nodes_txt.substr(pos, (next == std::string::npos) ? std::string::npos : next - pos);
            if (!tok.empty())
            {
                char dir = tok.back();
                char sym = (dir == '+') ? '>' : '<';
                oss << sym << tok.substr(0, tok.size() - 1);
            }
            if (next == std::string::npos)
                break;
            pos = next + 1;
        }

        return oss.str();
    }

    /**
     * @brief Convert the Walk to bGFA binary format
     * @param {ofstream} &output_file - Output binary file stream
     */
    void toBgfaBin(std::ofstream &output_file) const
    {
        if (!output_file.is_open())
        {
            throw std::runtime_error(
                "[Write Walk to bGFA File] Output file is not open");
        }

        // Store sample and haplotype
        output_file.write(reinterpret_cast<const char *>(&sample_id_), sizeof(uword_t));
        output_file.write(reinterpret_cast<const char *>(&haplotype_id_), sizeof(uword_t));
        output_file.write(reinterpret_cast<const char *>(&sequence_id_), sizeof(uword_t));

        path_.toBgfaBin(output_file);
    }

    /**
     * @brief Output stream operator for Walk
     * @param {ostream} &os - Output stream
     * @param {Walk} &walk - Walk to output
     * @return {ostream&} Reference to the output stream
     */
    friend std::ostream &operator<<(std::ostream &os, const Walk &walk)
    {
        os << "W";
        os << ", Sample: " << walk.sample_id_
           << ", Haplotype: " << walk.haplotype_id_
           << ", Sequence: " << walk.sequence_id_;

        const auto &nodes = walk.path_.getAllNodes();
        if (!nodes.empty())
        {
            os << "\n\t" << walk.dir_symbol(nodes[0]) << (nodes[0] >> 1);
            for (size_t i = 2; i < nodes.size(); ++i)
            {
                os << walk.dir_symbol(nodes[i]) << (nodes[i] >> 1);
            }
        }
        return os;
    }

private:
    char dir_symbol(uword_t node) const
    {
        return (node & 0x01) ? '>' : '<';
    };

    /**
     * @brief Compute storage size in uword_t units
     * @details Includes sample/haplotype fields plus path footprint
     */
    void setStoreSize()
    {
        // sample_id_, haplotype_id_ and sequence_id_
        store_size_ = 3;
        // plus path storage footprint
        store_size_ += path_.getStoreSize();
    }
};

// TODO:按照现在的设计，似乎需要一个判断segment的id不能超出 WORD_BIT>>2(?) 的机制
/**
 * @brief: GFA class representing the entire graph
 */
class GFA
{
private:
    uint64_t bp_num_; ///< Total base pairs

    uword_t seg_num_;  ///< Number of nodes
    uword_t link_num_; ///< Number of links
    uword_t path_num_; ///< Number of paths
    uword_t walk_num_; ///< Number of walks

    uword_t seg_size_;  ///< Size of segments
    uword_t link_size_; ///< Size of links
    uword_t path_size_; ///< Size of paths
    uword_t walk_size_; ///< Size of walks

    std::vector<Segment> segments_; ///< List of segments
    LinksSet links_;                ///< Set of links
    std::vector<Path> paths_;       ///< List of paths
    std::vector<Walk> walks_;       ///< List of walks
    NodeMap name_to_id_;            ///< Mapping from segment names to IDs
    NodeMap sn_name_to_id_;         ///< Mapping from SN names to numeric ids
    NodeMap w_sample_to_id_;        ///< Mapping from sample names to IDs
    NodeMap w_sequence_to_id_;      ///< Mapping from sample names to IDs

    // Section tag support for convert split mode
    uword_t pre_info_ = 0;
    // PATH_MODE 0x01=>inc、0x00=>direct
    bool store_sosr_ = true;
    bool sosr_compact_ = false;

public:
    /**
     * @brief Construct a new GFA object
     */
    GFA() : seg_num_(0), bp_num_(0), link_num_(0), path_num_(0), walk_num_(0),
            seg_size_(0), link_size_(0), path_size_(0), walk_size_(0)
    {
    }

    /**
     * @brief Loads GFA file (supports .gfa and .bgfa formats)
     * @param {string} &filename - Path to the GFA file
     */
    explicit GFA(const std::string &filename) : GFA()
    {
        // 获取文件后缀名
        size_t dot_pos = filename.find_last_of(".");
        if (dot_pos == std::string::npos || dot_pos == filename.length() - 1)
        {
            throw std::runtime_error(
                "[Read GFA File] Filename has no valid extension: " + filename);
        }
        std::string extension = filename.substr(dot_pos + 1);

        if (extension == "bgfa")
        {
            loadFromFileBgfa(filename);
        }
        else if (extension == "gfa")
        {
            loadFromFileGfa(filename);
        }
        else
        {
            throw std::logic_error(
                "[Read GFA File] Unsupported file format: ." + extension +
                ". Expected .gfa or .bgfa");
        }
    }

    /**
     * @brief Loads GFA file (supports .gfa and .bgfa formats)
     * @param {string} &filename - Path to the GFA file
     */
    explicit GFA(const std::string &filename, std::string path_mode) : GFA()
    {
        // 获取文件后缀名
        size_t dot_pos = filename.find_last_of(".");
        if (dot_pos == std::string::npos || dot_pos == filename.length() - 1)
        {
            throw std::runtime_error(
                "[Read GFA File] Filename has no valid extension: " + filename);
        }
        std::string extension = filename.substr(dot_pos + 1);

        if (path_mode == "inc")
        {
            pre_info_ |= PRE_INFO_PATHMODE_INC;
        }
        else if (path_mode == "direct")
        {
            pre_info_ &= (~PRE_INFO_PATHMODE_INC);
        }
        else
        {
            throw std::logic_error(
                "[Read GFA File] Unsupported path mode: " + path_mode +
                ". Expected 'inc' or 'direct'");
        }

        if (extension == "bgfa")
        {
            loadFromFileBgfa(filename);
        }
        else if (extension == "gfa")
        {
            loadFromFileGfa(filename);
        }
        else
        {
            throw std::logic_error(
                "[Read GFA File] Unsupported file format: ." + extension +
                ". Expected .gfa or .bgfa");
        }
    }

    /**
     * @brief Get the total number of base pairs in the GFA
     * @return {uint64_t} Total base pairs
     */
    uword_t getBpNum() const { return bp_num_; }

    /**
     * @brief Get the number of nodes in the GFA
     * @return {uword_t} Number of nodes
     */
    uword_t getNodeNum() const { return seg_num_; }

    /**
     * @brief Get the number of links in the GFA
     * @return {uword_t} Number of links
     */
    uword_t getLinkNum() const { return link_num_; }

    /**
     * @brief Get the number of paths in the GFA
     * @return {uword_t} Number of paths
     */
    uword_t getPathNum() const { return path_num_; }

    /**
     * @brief Get the number of walks in the GFA
     * @return {uword_t} Number of walks
     */
    uword_t getWalkNum() const { return walk_num_; }

    void setSosrStorage(bool store_sosr, bool compact_mode)
    {
        store_sosr_ = store_sosr;
        sosr_compact_ = compact_mode;
    }

    /**
     * @brief Get the size of segments in storage (in uword_t units)
     * @return {uword_t} Size of segments
     */
    uword_t getSegmentSize() const { return seg_size_; }

    /**
     * @brief Get the size of links in storage (in uword_t units)
     * @return {uword_t} Size of links
     */
    uword_t getLinkSize() const { return link_size_; }

    /**
     * @brief Get the size of paths in storage (in uword_t units)
     * @return {uword_t} Size of paths
     */
    uword_t getPathSize() const { return path_size_; }

    /**
     * @brief Get the size of walks in storage (in uword_t units)
     * @return {uword_t} Size of walks
     */
    uword_t getWalkSize() const { return walk_size_; }

    /**
     * @brief Get the segment index by its name
     * @param {string} &id
     * @return {uword_t} Segment index
     */
    uword_t getSegmentIndex(const std::string &id) const
    {
        auto it = name_to_id_.find(id);
        if (it == name_to_id_.end())
        {
            throw std::logic_error(
                "[Get Segment Index] Segment not found: " + id);
        }
        return it->second;
    }

    /**
     * @brief Get the segment by its index
     * @details Throws an out_of_range exception if the index is invalid
     *          It returns a constant reference to the Segment object
     * @param {uword_t} index - Segment index
     * @return {const Segment&} Segment object
     */
    const Segment &getSegmentByIndex(uword_t index) const
    {
        if (index >= segments_.size())
        {
            throw std::out_of_range(
                "[Get Segment By Index] Segment index out of range");
        }
        return segments_[index];
    }

    /**
     * @brief Get all segments in the GFA
     * @details It returns a constant reference to the vector of Segment objects
     * @return {const vector<Segment>&} Vector of all segments
     */
    const std::vector<Segment> &getAllSegments() const { return segments_; }

    /**
     * @brief Find segments by distance within a specified percentage difference
     * @param {uword_t} dis - Target distance
     * @param {uword_t} diff_precent - Allowed percentage difference
     * @return {vector<Segment>} Vector of matching segments (returned by value)
     */
    std::vector<Segment> findSegmentsByDistance(uword_t dis, uword_t diff_precent)
    {
        std::vector<Segment> result_segments;
        if (diff_precent > 100)
        {
            // clamp to 100% to avoid underflow/overflow
            diff_precent = 100;
        }
        uword_t lower_bound = (dis > 0) ? (dis - (dis * diff_precent) / 100) : 0;
        uword_t upper_bound = dis + (dis * diff_precent) / 100;

        for (const auto &segment : segments_)
        {
            uword_t seg_dis = segment.getDis();
            if (seg_dis >= lower_bound && seg_dis <= upper_bound)
            {
                result_segments.push_back(segment);
            }
            else if (seg_dis > upper_bound)
            {
                break;
            }
        }

        return result_segments;
    }

    /**
     * @brief Get the LinksSet of the GFA
     * @details It returns a constant reference to the LinksSet object
     * @return {const LinksSet&} LinksSet object
     */
    const LinksSet &getLinksSet() const { return links_; }

    /**
     * @brief Get all paths in the GFA
     * @details It returns a constant reference to the vector of Path objects
     * @return {const vector<Path>&} Vector of all paths
     */
    const std::vector<Path> &getAllPaths() const { return paths_; }

    /**
     * @brief Get all walks in the GFA
     * @details It returns a constant reference to the vector of Walk objects
     * @return {const vector<Walk>&} Vector of all walks
     */
    const std::vector<Walk> &getAllWalks() const { return walks_; }

    void briefStat() const
    {
        std::cout << "TotalBases: " << bp_num_ << std::endl;
        std::cout << "Nodes: " << seg_num_ << std::endl;
        std::cout << "Edges: " << link_num_ << std::endl;
        std::cout << "Paths: " << path_num_ << std::endl;
        std::cout << "Walks: " << walk_num_ << std::endl;
    }

    /**
     * @brief Print graph statistics:
     *  - Total bases
     *  - Node count and length distribution (<10, 10-30, 30-100, 100-300, 300-1000,
     *    1000-3000, 3000-10000, 10000-30000, >30000)
     *  - Edge count
     *  - Degree stats (avg/min/max/std) where degree = in-degree + out-degree
     *  - Degree distribution (count per degree value)
     */
    // TODO: Just for test, it is not optimized. It seems better to set some variables as class members to avoid recomputation
    void printStatistics() const
    {
        // Effective nodes (length>0)
        const size_t n_all = segments_.size();
        std::vector<uword_t> lengths;
        lengths.reserve(n_all);
        for (const auto &s : segments_)
        {
            if (s.getLength() > 0)
                lengths.push_back(s.getLength());
        }
        const size_t n_nodes = lengths.size();

        // Length bins
        // bins: <10, 10-30, 30-100, 100-300, 300-1000, 1000-3000, 3000-10000, 10000-30000, >30000
        std::array<uword_t, 9> len_bins{};
        len_bins.fill(0);
        auto bin_index = [](uword_t L) -> int
        {
            if (L < 10)
                return 0;
            if (L <= 30)
                return 1;
            if (L <= 100)
                return 2;
            if (L <= 300)
                return 3;
            if (L <= 1000)
                return 4;
            if (L <= 3000)
                return 5;
            if (L <= 10000)
                return 6;
            if (L <= 30000)
                return 7;
            return 8;
        };
        for (auto L : lengths)
        {
            len_bins[bin_index(L)]++;
        }

        // Edges
        const auto &all_links = links_.getAllLinks();
        uword_t edges = links_.getLinkNum();

        // Degrees per node (in + out) considering only effective nodes
        std::vector<uword_t> indeg(all_links.size(), 0), outdeg(all_links.size(), 0);
        for (size_t i = 0; i < all_links.size(); ++i)
        {
            outdeg[i] = static_cast<uword_t>(all_links[i].size());
            for (const auto &enc : all_links[i])
            {
                uword_t to = (enc >> LINK_DIR_BIT);
                if (to < indeg.size())
                    indeg[to]++;
            }
        }

        std::vector<double> degrees;
        degrees.reserve(n_all);
        for (size_t i = 0; i < segments_.size(); ++i)
        {
            if (segments_[i].getLength() == 0)
                continue;
            uword_t d = 0;
            if (i < indeg.size())
                d += indeg[i];
            if (i < outdeg.size())
                d += outdeg[i];
            degrees.push_back(static_cast<double>(d));
        }

        double avg = 0.0, mn = 0.0, mx = 0.0, stdv = 0.0;
        std::map<uword_t, uword_t> degree_hist;
        if (!degrees.empty())
        {
            double sum = std::accumulate(degrees.begin(), degrees.end(), 0.0);
            avg = sum / degrees.size();
            mn = degrees[0];
            mx = degrees[0];
            for (double d : degrees)
            {
                if (d < mn)
                    mn = d;
                if (d > mx)
                    mx = d;
                degree_hist[static_cast<uword_t>(d)]++;
            }
            double var_sum = 0.0;
            for (double d : degrees)
            {
                double diff = d - avg;
                var_sum += diff * diff;
            }
            stdv = std::sqrt(var_sum / degrees.size());
        }

        // Output
        std::cout << "[Stats] TotalBases: " << bp_num_ << std::endl;
        std::cout << "[Stats] Nodes: " << n_nodes << std::endl;
        std::cout << "[Stats] LengthBins:"
                  << "         <10       10-30      30-100     100-300"
                  << "    300-1000   1000-3000  3000-10000 10000-30000"
                  << "      >30000\n                   "
                  << std::setw(11) << len_bins[0] << ","
                  << std::setw(11) << len_bins[1] << ","
                  << std::setw(11) << len_bins[2] << ","
                  << std::setw(11) << len_bins[3] << ","
                  << std::setw(11) << len_bins[4] << ","
                  << std::setw(11) << len_bins[5] << ","
                  << std::setw(11) << len_bins[6] << ","
                  << std::setw(11) << len_bins[7] << ","
                  << std::setw(11) << len_bins[8] << std::endl;
        std::cout << "[Stats] Edges: " << edges << std::endl;
        std::cout << "[Stats] Degree avg: " << avg << ",\tmin: " << mn << ",\tmax: " << mx << ",\tstd: " << stdv << std::endl;
        std::cout << "[Stats] DegreeDistribution (deg=count): ";
        bool first = true;
        for (const auto &kv : degree_hist)
        {
            if (!first)
                std::cout << ", ";
            first = false;
            std::cout << kv.first << "=" << kv.second;
        }
        std::cout << std::endl;
    }

    /**
     * @brief Add a new segment to the GFA
     * @param {Segment} segment - Segment object to add
     * @return {uword_t} Index of the added segment
     */
    uword_t addSegment(Segment segment)
    {
        segment.setId(seg_num_);
        segment.setDis(0);
        segments_.push_back(segment);
        name_to_id_[std::to_string(segment.getId())] = seg_num_;
        links_.addNode();

        seg_num_++;
        bp_num_ += segment.getLength();
        seg_size_ += segment.getStoreSize();

        link_num_ = links_.getLinkNum();
        link_size_ = links_.getStoreSize();

        return seg_num_ - 1;
    }

    /**
     * @brief Adds a link using source and target node indices with directions
     * @param {uword_t} from_id - Source node index
     * @param {uword_t} to_id - Target node index
     * @param {string} from_dir - Direction from source node ('+' or '-')
     * @param {string} to_dir - Direction to target node ('+' or '-')
     */
    void addLink(uword_t from_id, uword_t to_id,
                 std::string from_dir, std::string to_dir)
    {
        links_.addLink(from_id, to_id, from_dir, to_dir);

        link_num_ = links_.getLinkNum();
        link_size_ = links_.getStoreSize();
    }

    /**
     * @brief Adds a link using source node index and encoded target value
     * @details This method adds a link using an encoded to_value
     * @param {uword_t} from_id - Source node index
     * @param {uword_t} to_value - Encoded target node index and direction
     */
    void addLink(uword_t from_id, uword_t to_value)
    {
        links_.addLink(from_id, to_value);

        link_num_ = links_.getLinkNum();
        link_size_ = links_.getStoreSize();
    }

    /**
     * @brief Add a new path to the GFA
     * @param {Path} path - Path object to add
     */
    void addPath(Path path)
    {
        paths_.push_back(path);

        path_num_++;
        path_size_ += path.getStoreSize();
    }

    /**
     * @brief Add a new walk to the GFA
     * @param {Walk} walk - Walk object to add
     */
    void addWalk(Walk walk)
    {
        walks_.push_back(walk);

        walk_num_++;
        walk_size_ += walk.getStoreSize();
    }

    /**
     * @brief Redirect all links pointing to old_to so they point to new_to (in-place).
     */
    void redirectLinksTo(uword_t old_to, uword_t new_to)
    {
        links_.redirectIncoming(old_to, new_to);

        link_num_ = links_.getLinkNum();
        link_size_ = links_.getStoreSize();
    }

    /**
     * @brief Move all outgoing links from old_from to new_from (in-place).
     */
    void redirectLinksFrom(uword_t old_from, uword_t new_from)
    {
        links_.redirectOutgoing(old_from, new_from);

        link_num_ = links_.getLinkNum();
        link_size_ = links_.getStoreSize();
    }

    /**
     * @brief Remove all outgoing/incoming links of a node.
     */
    void removeLinksOf(uword_t node_id)
    {
        links_.removeOutgoing(node_id);
        links_.removeIncoming(node_id);

        link_num_ = links_.getLinkNum();
        link_size_ = links_.getStoreSize();
    }

    /**
     * @brief Safely mark a segment as deleted without renumbering ids.
     *        The segment's sequence becomes empty and length is set to 0.
     */
    void markSegmentDeleted(uword_t node_id)
    {
        if (node_id >= segments_.size())
            throw std::out_of_range("[Mark Segment Deleted] Segment index out of range");

        bp_num_ -= segments_[node_id].getLength();
        seg_size_ -= segments_[node_id].getStoreSize();

        segments_[node_id].setSequence(std::vector<uword_t>());
        segments_[node_id].setDis(0);

        seg_size_ += segments_[node_id].getStoreSize();

        link_num_ = links_.getLinkNum();
        link_size_ = links_.getStoreSize();
    }

    /**
     * @brief Compact and renumber segments to sequential IDs, rebuilding links accordingly.
     *        Deleted/empty segments (length==0) are skipped.
     *        After renumbering, `name_to_id_` and `links_` are updated to match new IDs.
     */
    void compactAndRenumber()
    {
        // Build old->new id mapping, skipping empty segments
        std::unordered_map<uword_t, uword_t> id_map;
        std::vector<Segment> new_segments;
        new_segments.reserve(segments_.size());
        name_to_id_.clear();

        for (uword_t old_id = 0; old_id < segments_.size(); ++old_id)
        {
            const Segment &s = segments_[old_id];
            if (s.getLength() == 0)
                continue; // skip deleted
            Segment copy(0, s.getSequenceAsString());
            uword_t new_id = new_segments.size();
            copy.setId(new_id);
            copy.setDis(s.getDis());
            new_segments.push_back(copy);
            id_map[old_id] = new_id;
            name_to_id_[std::to_string(new_id)] = new_id;
        }

        // Rebuild links with new IDs
        LinksSet new_links;
        for (uword_t i = 0; i < new_segments.size(); ++i)
            new_links.addNode();

        const auto &old_links = links_.getAllLinks();
        for (uword_t from = 0; from < old_links.size(); ++from)
        {
            auto from_it = id_map.find(from);
            if (from_it == id_map.end())
                continue; // from deleted
            uword_t new_from = from_it->second;
            for (const auto &to_value : old_links[from])
            {
                uword_t old_to = (to_value >> LINK_DIR_BIT);
                uword_t dir_bits = (to_value & 0x03);
                auto to_it = id_map.find(old_to);
                if (to_it == id_map.end())
                {
                    // link to deleted node: drop
                    continue;
                }
                uword_t new_to = to_it->second;
                new_links.addLink(new_from, (new_to << LINK_DIR_BIT) | dir_bits);
            }
        }

        // Commit
        segments_.swap(new_segments);
        links_ = new_links;
        seg_num_ = segments_.size();
    }

    /**
     * @brief Merge linear unique pairs: if A has exactly one successor B,
     *        and B has exactly one predecessor A, merge A and B into one node.
     *        The merged node keeps the smaller id, and L/P/W are updated.
     * @return {uword_t} Number of successful merges
     */
    uword_t mergeLinearUniqueSuccessorPredecessor()
    {
        uword_t merged_count = 0;

        while (true)
        {
            const auto &all_links = links_.getAllLinks();
            if (all_links.empty())
                break;

            std::vector<uword_t> indeg(all_links.size(), 0);
            std::vector<uword_t> pred(all_links.size(), uword_t(-1));
            for (uword_t from = 0; from < all_links.size(); ++from)
            {
                for (const auto &enc : all_links[from])
                {
                    uword_t to = (enc >> LINK_DIR_BIT);
                    if (to < indeg.size())
                    {
                        indeg[to]++;
                        if (indeg[to] == 1)
                            pred[to] = from;
                    }
                }
            }

            bool found = false;
            uword_t a = 0, b = 0;

            for (uword_t from = 0; from < all_links.size(); ++from)
            {
                if (from >= segments_.size() || segments_[from].getLength() == 0)
                    continue;
                if (all_links[from].size() != 1)
                    continue;

                uword_t enc = *(all_links[from].begin());
                uword_t to = (enc >> LINK_DIR_BIT);
                if (to >= segments_.size() || segments_[to].getLength() == 0)
                    continue;
                if (to == from)
                    continue;

                if (to < indeg.size() && indeg[to] == 1 && pred[to] == from)
                {
                    a = from;
                    b = to;
                    found = true;
                    break;
                }
            }

            if (!found)
                break;

            uword_t keep = std::min(a, b);
            uword_t drop = std::max(a, b);

            std::string merged_seq = segments_[a].getSequenceAsString() +
                                     segments_[b].getSequenceAsString();

            segments_[keep].setSequenceFromString(merged_seq);
            segments_[keep].setDis(std::min(segments_[a].getDis(), segments_[b].getDis()));

            segments_[drop].setSequence(std::vector<uword_t>());
            segments_[drop].setDis(0);

            auto old_links = links_.getAllLinks();
            LinksSet new_links;
            for (uword_t i = 0; i < segments_.size(); ++i)
                new_links.addNode();

            for (uword_t from = 0; from < old_links.size(); ++from)
            {
                for (const auto &enc : old_links[from])
                {
                    uword_t to = (enc >> LINK_DIR_BIT);
                    uword_t dir = (enc & 0x03);

                    if (from == a && to == b)
                        continue;

                    uword_t new_from = (from == drop) ? keep : from;
                    uword_t new_to = (to == drop) ? keep : to;

                    if (new_from >= segments_.size() || new_to >= segments_.size())
                        continue;

                    new_links.addLink(new_from, (new_to << LINK_DIR_BIT) | dir);
                }
            }
            links_ = new_links;

            std::vector<Path> new_paths;
            new_paths.reserve(paths_.size());
            for (const auto &path : paths_)
            {
                std::string line = path.toGfaTxt();
                auto nodes = parsePathNodesFromLine(line);
                auto nodes_new = remapMergedNodes(nodes, a, b, keep, drop);
                if (nodes_new.empty())
                    continue;

                std::string body = buildPathBody(nodes_new);
                PathMode mode = path.getMode();
                new_paths.emplace_back(Path(path.getId(), body, name_to_id_, mode));
            }
            paths_.swap(new_paths);

            std::vector<Walk> new_walks;
            new_walks.reserve(walks_.size());
            for (const auto &walk : walks_)
            {
                const Path &wpath = walk.getPath();
                std::string line = wpath.toGfaTxt();
                auto nodes = parsePathNodesFromLine(line);
                auto nodes_new = remapMergedNodes(nodes, a, b, keep, drop);
                if (nodes_new.empty())
                    continue;

                std::string body = buildPathBody(nodes_new);
                PathMode mode = wpath.getMode();
                Path new_path(0, body, name_to_id_, mode);
                new_walks.emplace_back(Walk(walk.getSampleId(), walk.getHaplotypeId(),
                                            walk.getSequenceId(), new_path));
            }
            walks_.swap(new_walks);

            refreshMetaStats();
            merged_count++;
        }

        return merged_count;
    }

    uword_t mergeUniqueSuccessorByReverseTable()
    {
        const uword_t n = static_cast<uword_t>(segments_.size());
        if (n == 0)
            return 0;

        const auto &all_links = links_.getAllLinks();
        if (all_links.empty())
            return 0;

        std::vector<uint8_t> valid(n, 0);
        for (uword_t i = 0; i < n; ++i)
        {
            valid[i] = (segments_[i].getLength() > 0) ? 1 : 0;
        }

        std::vector<uword_t> indeg(n, 0), outdeg(n, 0);
        for (uword_t from = 0; from < all_links.size(); ++from)
        {
            if (from >= n || !valid[from])
                continue;

            outdeg[from] = static_cast<uword_t>(all_links[from].size());
            for (const auto &enc : all_links[from])
            {
                const uword_t to = (enc >> LINK_DIR_BIT);
                if (to < n && valid[to])
                    indeg[to]++;
            }
        }

        // Successor table for merge-eligible edges only:
        // A+ -> B+ and outdeg[A]==1 and indeg[B]==1.
        std::vector<int64_t> succ(n, -1);
        for (uword_t from = 0; from < all_links.size(); ++from)
        {
            if (from >= n || !valid[from])
                continue;
            if (outdeg[from] != 1)
                continue;

            const uword_t enc = *(all_links[from].begin());
            const uword_t to = (enc >> LINK_DIR_BIT);
            const uword_t dir = (enc & 0x03);
            if (dir != 0x03)
                continue;
            if (to >= n || !valid[to] || to == from)
                continue;
            if (indeg[to] != 1)
                continue;

            succ[from] = static_cast<int64_t>(to);
        }

        // indegree in successor-subgraph (outdegree <=1 by construction)
        std::vector<uword_t> indeg_s(n, 0);
        for (uword_t u = 0; u < n; ++u)
        {
            if (succ[u] >= 0)
                indeg_s[static_cast<uword_t>(succ[u])]++;
        }

        std::vector<uint8_t> visited(n, 0);
        std::vector<std::vector<uword_t>> chains;

        auto collect_linear_from = [&](uword_t start)
        {
            std::vector<uword_t> chain;
            uword_t x = start;
            while (x < n && !visited[x] && valid[x])
            {
                chain.push_back(x);
                visited[x] = 1;
                if (succ[x] < 0)
                    break;
                x = static_cast<uword_t>(succ[x]);
            }
            if (chain.size() >= 2)
                chains.push_back(std::move(chain));
        };

        // Phase 1: start from heads in successor-subgraph
        for (uword_t u = 0; u < n; ++u)
        {
            if (!valid[u] || visited[u] || succ[u] < 0)
                continue;
            if (indeg_s[u] == 1)
                continue;
            collect_linear_from(u);
        }

        // Phase 2: remaining unvisited nodes belong to pure cycles
        for (uword_t u = 0; u < n; ++u)
        {
            if (!valid[u] || visited[u] || succ[u] < 0)
                continue;

            std::vector<uword_t> cyc;
            std::unordered_map<uword_t, size_t> pos;
            uword_t x = u;
            while (x < n && !visited[x] && valid[x])
            {
                if (pos.find(x) != pos.end())
                    break;
                pos[x] = cyc.size();
                cyc.push_back(x);
                if (succ[x] < 0)
                    break;
                x = static_cast<uword_t>(succ[x]);
            }

            for (auto node : cyc)
                visited[node] = 1;

            auto it = pos.find(x);
            if (it != pos.end())
            {
                std::vector<uword_t> pure_cycle(cyc.begin() + static_cast<ptrdiff_t>(it->second), cyc.end());
                if (pure_cycle.size() >= 2)
                    chains.push_back(std::move(pure_cycle));
            }
            else if (cyc.size() >= 2)
            {
                chains.push_back(std::move(cyc));
            }
        }

        if (chains.empty())
            return 0;

        uword_t merged_count = 0;
        for (const auto &chain : chains)
        {
            if (chain.size() >= 2)
                merged_count += static_cast<uword_t>(chain.size() - 1);
        }

        // old id -> new id, assigned in rebuild phase
        std::vector<uword_t> old2new(n, uword_t(-1));
        std::vector<Segment> new_segments;
        new_segments.reserve(n);

        // Merge each chain into one node (chain order defines sequence order)
        for (const auto &chain : chains)
        {
            std::string merged_seq;
            uword_t min_dis = std::numeric_limits<uword_t>::max();
            for (auto node : chain)
            {
                merged_seq += segments_[node].getSequenceAsString();
                min_dis = std::min(min_dis, segments_[node].getDis());
            }

            uword_t new_id = static_cast<uword_t>(new_segments.size());
            Segment merged_seg(new_id, merged_seq);
            merged_seg.setDis(min_dis == std::numeric_limits<uword_t>::max() ? 0 : min_dis);
            new_segments.push_back(merged_seg);

            for (auto node : chain)
                old2new[node] = new_id;
        }

        // Keep all remaining valid nodes as singletons
        for (uword_t old_id = 0; old_id < n; ++old_id)
        {
            if (!valid[old_id])
                continue;
            if (old2new[old_id] != uword_t(-1))
                continue;

            uword_t new_id = static_cast<uword_t>(new_segments.size());
            Segment copy_seg(new_id, segments_[old_id].getSequenceAsString());
            copy_seg.setDis(segments_[old_id].getDis());
            new_segments.push_back(copy_seg);
            old2new[old_id] = new_id;
        }

        // Rebuild name map
        name_to_id_.clear();
        for (uword_t i = 0; i < new_segments.size(); ++i)
            name_to_id_[std::to_string(i)] = i;

        // Rebuild links in one pass
        LinksSet new_links;
        for (uword_t i = 0; i < new_segments.size(); ++i)
            new_links.addNode();

        for (uword_t from = 0; from < all_links.size(); ++from)
        {
            if (from >= n || !valid[from])
                continue;
            const uword_t new_from = old2new[from];
            if (new_from == uword_t(-1))
                continue;

            for (const auto &enc : all_links[from])
            {
                const uword_t to = (enc >> LINK_DIR_BIT);
                const uword_t dir = (enc & 0x03);
                if (to >= n || !valid[to])
                    continue;
                const uword_t new_to = old2new[to];
                if (new_to == uword_t(-1))
                    continue;
                if (new_from == new_to)
                    continue; // collapsed inside merged node

                new_links.addLink(new_from, (new_to << LINK_DIR_BIT) | dir);
            }
        }

        // Rebuild paths and walks once, with remap + normalization
        std::vector<Path> new_paths;
        new_paths.reserve(paths_.size());
        for (const auto &path : paths_)
        {
            const std::string line = path.toGfaTxt();
            const auto nodes = parsePathNodesFromLine(line);
            const auto mapped = remapNodesByMapNormalized(nodes, old2new);
            if (mapped.empty())
                continue;

            const std::string body = buildPathBody(mapped);
            const PathMode mode = path.getMode();
            new_paths.emplace_back(Path(path.getId(), body, name_to_id_, mode));
        }

        std::vector<Walk> new_walks;
        new_walks.reserve(walks_.size());
        for (const auto &walk : walks_)
        {
            const Path &wpath = walk.getPath();
            const std::string line = wpath.toGfaTxt();
            const auto nodes = parsePathNodesFromLine(line);
            const auto mapped = remapNodesByMapNormalized(nodes, old2new);
            if (mapped.empty())
                continue;

            const std::string body = buildPathBody(mapped);
            const PathMode mode = wpath.getMode();
            Path new_path(0, body, name_to_id_, mode);
            new_walks.emplace_back(Walk(walk.getSampleId(), walk.getHaplotypeId(),
                                        walk.getSequenceId(), new_path));
        }

        // Commit
        segments_.swap(new_segments);
        links_ = new_links;
        paths_.swap(new_paths);
        walks_.swap(new_walks);
        refreshMetaStats();

        return merged_count;
    }

    /**
     * @brief  Set the distances for all segments in the GFA
     */
    void setSegmentsDis()
    {
        std::queue<uword_t> segements_queue;
        std::unordered_set<uword_t> has_predecessor;

        for (const auto &link : links_.getAllLinks())
        {
            for (const auto &to_link : link)
            {
                has_predecessor.insert(to_link >> LINK_DIR_BIT);
            }
        }
        for (uword_t i = 0; i < seg_num_; ++i)
        {
            if (has_predecessor.find(i) == has_predecessor.end())
            {
                segements_queue.push(i);
                segments_[i].setDis(0);
            }
        }

        while (!segements_queue.empty())
        {
            uword_t current_id = segements_queue.front();
            segements_queue.pop();
            uword_t current_dis = segments_[current_id].getDis();
            for (const auto &to_link : links_.getLinksFromNode(current_id))
            {
                uword_t to_id = to_link >> LINK_DIR_BIT;
                uword_t dir = to_link & 0x03;

                if (segments_[to_id].getDis() == 0 ||
                    segments_[to_id].getDis() >
                        current_dis + segments_[current_id].getLength())
                {
                    segments_[to_id].setDis(
                        current_dis + segments_[current_id].getLength());
                    segements_queue.push(to_id);
                }
            }
        }

        // It can not be deleted cause used in graph merge.
        // std::sort(segments_.begin(), segments_.end());
    }

    /**
     * @brief Check if a segment is in the GFA graph
     * @details It has two overloads:
     * @param {Segment} seg - Segment to check
     * @return {uword_t} Segment ID if found, -1 otherwise
     */
    uword_t isSegmentInGraph(Segment seg)
    {
        for (const auto &s : segments_)
        {
            if (s == seg)
                return s.getId();
        }
        return uword_t(-1);
    }

    /**
     * @brief Output stream operator for GFA
     * @return {ostream&} Reference to the output stream
     */
    friend std::ostream &operator<<(std::ostream &os, const GFA &gfa)
    {
        for (const auto &segment : gfa.segments_)
        {
            os << segment << std::endl;
        }

        os << gfa.links_ << std::endl;

        for (const auto &path : gfa.paths_)
        {
            os << path << std::endl;
        }

        for (const auto &walk : gfa.walks_)
        {
            os << walk << std::endl;
        }

        return os;
    }

    /**
     * @brief Write the GFA graph to a GFA format file
     * @param {string} &filename - Output filename
     */
    void write2GFA(const std::string &filename) const
    {
        std::ofstream file(filename);
        if (!file.is_open())
        {
            throw std::runtime_error(
                "[Write to GFA File] Failed to open file for writing: " +
                filename);
        }

        file << "H\tVN:Z:1.0" << std::endl;

        // Write segments
        for (const auto &segment : segments_)
        {
            if (segment.getLength() == 0)
                continue; // skip deleted segments
            file << segment.toGfaTxt() << std::endl;
        }

        file << links_.toGfaTxt();

        for (const auto &path : paths_)
        {
            file << path.toGfaTxt() << std::endl;
        }

        for (const auto &walk : walks_)
        {
            file << walk.toGfaTxt() << std::endl;
        }

        file.close();
    }

    /**
     * @brief Write the GFA graph to a bGFA binary format file
     * @param {string} &filename - Output filename
     */
    void write2Bgfa(const std::string &filename, bool store_segment_id = true)
    {
        std::ofstream ofile(filename, std::ios::binary);
        if (!ofile.is_open())
        {
            throw std::runtime_error(
                "[Write to bGFA File] Failed to open file for writing: " +
                filename);
        }

        // setPreInfo(ofile);

        pre_info_ &= PRE_INFO_PATHMODE_INC;
        if (!store_segment_id)
            pre_info_ |= PRE_INFO_SEGMENT_NO_ID;
        else
            pre_info_ &= (~PRE_INFO_SEGMENT_NO_ID);
        if (store_sosr_)
        {
            pre_info_ |= PRE_INFO_HAS_SOSR;
            if (!sosr_compact_)
                pre_info_ |= PRE_INFO_SOSR_UNCOMPRESSED;
            else
                pre_info_ &= (~PRE_INFO_SOSR_UNCOMPRESSED);
        }
        else
        {
            pre_info_ &= (~PRE_INFO_HAS_SOSR);
            pre_info_ &= (~PRE_INFO_SOSR_UNCOMPRESSED);
        }
        seg_size_ = 0;
        link_size_ = 0;
        path_size_ = 0;
        walk_size_ = 0;

        ofile.write(reinterpret_cast<const char *>(&pre_info_),
                    sizeof(uword_t));
        ofile.write(reinterpret_cast<const char *>(&seg_size_),
                    sizeof(uword_t));
        ofile.write(reinterpret_cast<const char *>(&link_size_),
                    sizeof(uword_t));
        ofile.write(reinterpret_cast<const char *>(&path_size_),
                    sizeof(uword_t));
        ofile.write(reinterpret_cast<const char *>(&walk_size_),
                    sizeof(uword_t));

        // Write segments
        for (auto &segment : segments_)
        {
            segment.setSosrStorage(store_sosr_, sosr_compact_);
            segment.toBgfaBin(ofile, store_segment_id);
            seg_size_ += segment.getStoreSize(store_segment_id);
        }

        {
            links_.toBgfaBin(ofile);
            link_size_ += links_.getStoreSize();
        }

        for (auto &path : paths_)
        {
            path.toBgfaBin(ofile);
            path_size_ += path.getStoreSize();
        }

        for (auto &walk : walks_)
        {
            walk.toBgfaBin(ofile);
            walk_size_ += walk.getStoreSize();
        }

        ofile.close();

        std::fstream file(filename,
                          std::ios::binary | std::ios::in | std::ios::out);
        file.seekp((WORD_BIT / 8), std::ios::beg);
        file.write(reinterpret_cast<const char *>(&seg_size_), sizeof(uword_t));
        file.write(reinterpret_cast<const char *>(&link_size_), sizeof(uword_t));
        file.write(reinterpret_cast<const char *>(&path_size_), sizeof(uword_t));
        file.write(reinterpret_cast<const char *>(&walk_size_), sizeof(uword_t));
        // If section tags are enabled, write four uword tags after sizes

        file.flush();
        file.close();
    }

    /**
     * @brief Write split binary files for S/L/P/W sections and a main BGFA with tags.
     * @param path_noext Base output path without extension
     */
    void write2BgfaSplit(const std::string &path_noext, bool store_segment_id = true)
    {
        const std::string s_file = path_noext + ".bgfas";
        const std::string l_file = path_noext + ".bgfal";
        const std::string p_file = path_noext + ".bgfap";
        const std::string w_file = path_noext + ".bgfaw";
        const std::string main_bgfa = path_noext + ".bgfa";

        // setPreInfo(ofile);
        // Generate random uword tags
        std::random_device rd;
        std::mt19937_64 gen(rd());
        using limit_uword = std::numeric_limits<uword_t>;
        std::uniform_int_distribution<uint64_t> dist(0ull, static_cast<uint64_t>(limit_uword::max()));
        uword_t s_tag = static_cast<uword_t>(dist(gen));
        uword_t l_tag = static_cast<uword_t>(dist(gen));
        uword_t p_tag = static_cast<uword_t>(dist(gen));
        uword_t w_tag = static_cast<uword_t>(dist(gen));

        pre_info_ = PRE_INFO_PATHMODE_INC;
        if (!store_segment_id)
            pre_info_ |= PRE_INFO_SEGMENT_NO_ID;
        if (store_sosr_)
        {
            pre_info_ |= PRE_INFO_HAS_SOSR;
            if (!sosr_compact_)
                pre_info_ |= PRE_INFO_SOSR_UNCOMPRESSED;
            else
                pre_info_ &= (~PRE_INFO_SOSR_UNCOMPRESSED);
        }
        else
        {
            pre_info_ &= (~PRE_INFO_HAS_SOSR);
            pre_info_ &= (~PRE_INFO_SOSR_UNCOMPRESSED);
        }

        {
            std::ofstream ofile(main_bgfa, std::ios::binary);
            if (!ofile.is_open())
            {
                throw std::runtime_error(
                    "[Write to bGFA File] Failed to open file for writing: " +
                    main_bgfa);
            }

            ofile.write(reinterpret_cast<const char *>(&pre_info_),
                        sizeof(uword_t));
            ofile.write(reinterpret_cast<const char *>(&s_tag),
                        sizeof(uword_t));
            ofile.write(reinterpret_cast<const char *>(&l_tag),
                        sizeof(uword_t));
            ofile.write(reinterpret_cast<const char *>(&p_tag),
                        sizeof(uword_t));
            ofile.write(reinterpret_cast<const char *>(&w_tag),
                        sizeof(uword_t));
        }

        // Write S
        {
            std::ofstream ofs(s_file, std::ios::binary);
            if (!ofs.is_open())
                throw std::runtime_error(std::string("[Split] Failed to open: ") + s_file);
            ofs.write(reinterpret_cast<const char *>(&s_tag), sizeof(uword_t));
            for (auto &segment : segments_)
            {
                segment.setSosrStorage(store_sosr_, sosr_compact_);
                segment.toBgfaBin(ofs, store_segment_id);
            }
            ofs.close();
        }
        // Write L
        {
            std::ofstream ofs(l_file, std::ios::binary);
            if (!ofs.is_open())
                throw std::runtime_error(std::string("[Split] Failed to open: ") + l_file);
            ofs.write(reinterpret_cast<const char *>(&l_tag), sizeof(uword_t));
            links_.toBgfaBin(ofs);
            ofs.close();
        }
        // Write P
        {
            std::ofstream ofs(p_file, std::ios::binary);
            if (!ofs.is_open())
                throw std::runtime_error(std::string("[Split] Failed to open: ") + p_file);
            ofs.write(reinterpret_cast<const char *>(&p_tag), sizeof(uword_t));
            for (auto &path : paths_)
                path.toBgfaBin(ofs);
            ofs.close();
        }
        // Write W
        {
            std::ofstream ofs(w_file, std::ios::binary);
            if (!ofs.is_open())
                throw std::runtime_error(std::string("[Split] Failed to open: ") + w_file);
            ofs.write(reinterpret_cast<const char *>(&w_tag), sizeof(uword_t));
            for (auto &walk : walks_)
                walk.toBgfaBin(ofs);
            ofs.close();
        }
    }

private:
    std::vector<std::pair<uword_t, char>> remapNodesByMapNormalized(
        const std::vector<std::pair<uword_t, char>> &nodes,
        const std::vector<uword_t> &old2new) const
    {
        std::vector<std::pair<uword_t, char>> out;
        out.reserve(nodes.size());

        for (const auto &nd : nodes)
        {
            const uword_t old_id = nd.first;
            const char dir = nd.second;
            if (old_id >= old2new.size())
                continue;
            const uword_t new_id = old2new[old_id];
            if (new_id == uword_t(-1))
                continue;

            if (!out.empty() && out.back().first == new_id)
            {
                // collapse consecutive duplicate mapped nodes
                continue;
            }
            out.emplace_back(new_id, dir);
        }

        return out;
    }

    std::vector<std::pair<uword_t, char>> parsePathNodesFromLine(const std::string &line) const
    {
        std::vector<std::pair<uword_t, char>> nodes;

        size_t t1 = line.find('\t');
        if (t1 == std::string::npos)
            return nodes;
        size_t t2 = line.find('\t', t1 + 1);
        if (t2 == std::string::npos)
            return nodes;
        size_t t3 = line.find('\t', t2 + 1);

        std::string body = (t3 == std::string::npos)
                               ? line.substr(t2 + 1)
                               : line.substr(t2 + 1, t3 - (t2 + 1));

        std::istringstream iss(body);
        std::string token;
        while (std::getline(iss, token, ','))
        {
            if (token.size() < 2)
                continue;
            char dir = token.back();
            if (dir != '+' && dir != '-')
                continue;

            std::string id_str = token.substr(0, token.size() - 1);
            uword_t id = static_cast<uword_t>(std::stoull(id_str));
            nodes.emplace_back(id, dir);
        }

        return nodes;
    }

    std::vector<std::pair<uword_t, char>> remapMergedNodes(
        const std::vector<std::pair<uword_t, char>> &nodes,
        uword_t a, uword_t b, uword_t keep, uword_t drop) const
    {
        std::vector<std::pair<uword_t, char>> out;
        out.reserve(nodes.size());

        for (size_t i = 0; i < nodes.size(); ++i)
        {
            uword_t old_id = nodes[i].first;
            char dir = nodes[i].second;

            if (i > 0 && nodes[i - 1].first == a && old_id == b)
            {
                continue;
            }

            uword_t new_id = (old_id == drop) ? keep : old_id;
            out.emplace_back(new_id, dir);
        }

        return out;
    }

    std::string buildPathBody(const std::vector<std::pair<uword_t, char>> &nodes) const
    {
        std::ostringstream oss;
        for (size_t i = 0; i < nodes.size(); ++i)
        {
            if (i > 0)
                oss << ",";
            oss << nodes[i].first << nodes[i].second;
        }
        return oss.str();
    }

    void refreshMetaStats()
    {
        bp_num_ = 0;
        seg_size_ = 0;
        for (auto &segment : segments_)
        {
            if (segment.getLength() > 0)
                bp_num_ += segment.getLength();
            seg_size_ += segment.getStoreSize();
        }

        seg_num_ = segments_.size();
        link_num_ = links_.getLinkNum();
        link_size_ = links_.getStoreSize();

        path_num_ = paths_.size();
        path_size_ = 0;
        for (auto &path : paths_)
        {
            path_size_ += path.getStoreSize();
        }

        walk_num_ = walks_.size();
        walk_size_ = 0;
        for (auto &walk : walks_)
        {
            walk_size_ += walk.getStoreSize();
        }
    }

    /**
     * @brief Load GFA graph from a GFA format file
     * @param {string} &filename - Input filename
     */
    void loadFromFileGfa(const std::string &filename)
    {

        std::ifstream file(filename);
        if (!file.is_open())
        {
            throw std::runtime_error(
                "[Load From GFA File] Failed to open file: " + filename);
        }

        std::string line;
        while (std::getline(file, line))
        {
            if (line.empty() || line[0] == '#')
                continue;

            if (line[0] == 'S')
            {
                parseSegmentLine(line);
            }
            else if (line[0] == 'L')
            {
                parseLinkLine(line);
            }
            else if (line[0] == 'P')
            {
                parsePathLine(line);
            }
            else if (line[0] == 'W')
            {
                parseWalkLine(line);
            }
        }

        links_.resizeLinks(seg_num_);

#ifdef INFO_DEBUG
        std::cout << "[Read GFA] Finished loading GFA file: " << filename << std::endl;
        briefStat();
        std::cout << "LinkSetSize:" << links_.getAllLinks().size() << std::endl;
#endif
        link_size_ = links_.getStoreSize();

        // pre_info_ = 0x01;

        // // Ensure links row size matches node count for serialization
        // links_.resizeLinks(node_num_);

        file.close();
    }

    /**
     * @brief Load GFA graph from a bGFA binary format file
     * @param {string} &filename - Input filename
     */
    void loadFromFileBgfa(const std::string &filename)
    {
        std::ifstream file(filename, std::ios::binary);
        uword_t curr_pos = 0;

        if (!file.is_open())
        {
            throw std::runtime_error(
                "[Load From bGFA File] Failed to open file: " + filename);
        }

        file.read(reinterpret_cast<char *>(&pre_info_), sizeof(pre_info_));
        file.read(reinterpret_cast<char *>(&seg_size_), sizeof(seg_size_));
        file.read(reinterpret_cast<char *>(&link_size_), sizeof(link_size_));
        file.read(reinterpret_cast<char *>(&path_size_), sizeof(path_size_));
        file.read(reinterpret_cast<char *>(&walk_size_), sizeof(walk_size_));

        store_sosr_ = ((pre_info_ & PRE_INFO_HAS_SOSR) != 0);
        sosr_compact_ = ((pre_info_ & PRE_INFO_SOSR_UNCOMPRESSED) == 0);

        // curr_pos;

        while (true)
        {
            // Segment record
            if (curr_pos < seg_size_)
            {
                parseBinSegmentField(file, curr_pos, seg_size_);
            }
            // Link record
            else if (curr_pos < seg_size_ + link_size_)
            {
                parseBinLinkField(file, curr_pos, seg_size_ + link_size_);
            }
            // Path record
            else if (curr_pos < seg_size_ + link_size_ + path_size_)
            {
                parseBinPathField(file, curr_pos, seg_size_ + link_size_ + path_size_);
            }
            // Walk record
            else if (curr_pos < seg_size_ + link_size_ + path_size_ + walk_size_)
            {
                parseBinWalkField(file, curr_pos,
                                  seg_size_ + link_size_ + path_size_ + walk_size_);
            }
            else
                break;

            if (file.eof())
                break;
        }
    }

    /**
     * @brief Parse a segment line from GFA file
     * @param {string} &line - Line from GFA file
     */
    void parseSegmentLine(const std::string &line)
    {
        std::stringstream ss(line);
        std::string type, id, sequence;
        uword_t so = 0, sr = 0, sn_id = 0;
        bool has_sn = false;

        ss >> type >> id >> sequence;

        std::string token;
        while (ss >> token)
        {
            // SR:i:<rank>
            const std::string sr_prefix = "SR:i:";
            if (token.compare(0, sr_prefix.size(), sr_prefix) == 0)
            {
                const std::string val = token.substr(sr_prefix.size());
                // If invalid, ignore and keep default 0
                try
                {
                    sr = static_cast<uword_t>(std::stoull(val));
                }
                catch (...)
                {
                }
                continue;
            }

            // SO:i:<offset>
            const std::string so_prefix = "SO:i:";
            if (token.compare(0, so_prefix.size(), so_prefix) == 0)
            {
                const std::string val = token.substr(so_prefix.size());
                // If invalid, ignore and keep default 0
                try
                {
                    so = static_cast<uword_t>(std::stoull(val));
                }
                catch (...)
                {
                }
                continue;
            }
            // Ignore other optional fields (e.g., LN:i:<len>, etc.)
            const std::string sn_prefix = "SN:Z:";
            if (token.compare(0, sn_prefix.size(), sn_prefix) == 0)
            {
                const std::string name = token.substr(sn_prefix.size());
                if (!name.empty())
                {
                    auto it = sn_name_to_id_.find(name);
                    if (it == sn_name_to_id_.end())
                    {
                        sn_id = static_cast<uword_t>(sn_name_to_id_.size());
                        sn_name_to_id_[name] = sn_id;
                    }
                    else
                    {
                        sn_id = it->second;
                    }
                    has_sn = true;
                }
                continue;
            }
        }

        name_to_id_[id] = seg_num_;
        Segment segment(seg_num_, sequence, so, sr);
        if (has_sn)
            segment.setSnId(sn_id);

        seg_num_++;
        bp_num_ += segment.getLength();
        seg_size_ += segment.getStoreSize();
        segments_.push_back(segment);
    }

    /**
     * @brief Parse a binary segment record from bGFA file
     * @param {ifstream} &file - Input binary file stream
     * @param {uword_t} &curr_pos - Current position in the file (in uword_t units)
     * @param {uword_t} limit - Limit position for segment records
     */
    void parseBinSegmentField(std::ifstream &file, uword_t &curr_pos, const uword_t limit)
    {
        uword_t node_num = seg_num_, sosr_word = 0, sn_word = 0, l_str = 0;
        const bool no_segment_id = ((pre_info_ & PRE_INFO_SEGMENT_NO_ID) != 0);
        const bool has_sosr = ((pre_info_ & PRE_INFO_HAS_SOSR) != 0);
        const bool sosr_uncompressed = ((pre_info_ & PRE_INFO_SOSR_UNCOMPRESSED) != 0);
        if (!no_segment_id)
        {
            file.read(reinterpret_cast<char *>(&node_num), sizeof(node_num));
            curr_pos += 1;
        }
        if (has_sosr)
        {
            file.read(reinterpret_cast<char *>(&sosr_word), sizeof(sosr_word));
            curr_pos += 1;
            if (sosr_uncompressed)
            {
                file.read(reinterpret_cast<char *>(&sn_word), sizeof(sn_word));
                curr_pos += 1;
            }
        }
        file.read(reinterpret_cast<char *>(&l_str), sizeof(l_str));
        curr_pos += 1;

        if (!file)
        {
            throw std::runtime_error(
                "[Load From bGFA File] Error reading S record header.");
        }
        uword_t so = 0, sr = 0, sn_id = 0;
        if (has_sosr)
        {
            if (sosr_uncompressed)
            {
                sr = (sosr_word & 0x7F);
                so = (sosr_word >> 7) & 0x01FFFFFF;
                sn_id = (sn_word & 0xFFFFFFFF);
            }
            else
            {
                sr = (sosr_word & 0x07);
                so = (sosr_word >> 3) & 0x003FFFFF;
                sn_id = (sosr_word >> 25) & 0x7F;
            }
        }

        if (no_segment_id)
        {
            const uword_t tag2 = (l_str >> (WORD_BIT - 2)) & 0x3;

            if (tag2 == 0)
            {
                std::vector<uword_t> seq_data;
                uint8_t baseinElement = WORD_BIT / 2;
                size_t seq_num = (l_str + baseinElement - 1) / baseinElement;
                seq_data.resize(seq_num);
                if (seq_num > 0)
                {
                    file.read(reinterpret_cast<char *>(seq_data.data()),
                              seq_num * sizeof(uword_t));
                    if (!file)
                    {
                        throw std::runtime_error("[Load From bGFA File] Error reading sequence data in S record.");
                    }
                }
                curr_pos += seq_num;
                Segment segment(node_num, so, sr, l_str, seq_data);
                segment.setSosrStorage(has_sosr, !sosr_uncompressed);
                if (has_sosr)
                    segment.setSnId(sn_id);
                bp_num_ += segment.getLength();
                segments_.push_back(segment);
            }
            else if (tag2 == 2)
            {
                const uword_t len = (l_str >> (WORD_BIT - 6)) & 0x0F;
                const uword_t payload = l_str & ((uword_t(1) << (WORD_BIT - 6)) - 1);
                std::string seq;
                seq.reserve(len);
                for (uword_t i = 0; i < len; ++i)
                {
                    const uword_t shift = (WORD_BIT - 8) - 2 * i;
                    const uint8_t enc = static_cast<uint8_t>((payload >> shift) & 0x3);
                    seq.push_back(GlobalVariant::bin2char[enc]);
                }
                Segment segment(node_num, seq, so, sr);
                segment.setSosrStorage(has_sosr, !sosr_uncompressed);
                if (has_sosr)
                    segment.setSnId(sn_id);
                bp_num_ += segment.getLength();
                segments_.push_back(segment);
            }
            else if (tag2 == 3)
            {
                uword_t payload_lo = 0;
                file.read(reinterpret_cast<char *>(&payload_lo), sizeof(payload_lo));
                if (!file)
                {
                    throw std::runtime_error("[Load From bGFA File] Error reading compact-2 sequence payload in S record.");
                }
                curr_pos += 1;

                const uword_t len = (l_str >> (WORD_BIT - 8)) & 0x3F;
                const uword_t hi_bits = WORD_BIT - 8;
                const uword_t payload_hi = l_str & ((uword_t(1) << hi_bits) - 1);

                std::string seq;
                seq.reserve(len);
                for (uword_t i = 0; i < len; ++i)
                {
                    const uword_t p0 = 2 * i;
                    const uword_t p1 = p0 + 1;

                    const uint8_t b0 = (p0 < hi_bits)
                                           ? static_cast<uint8_t>((payload_hi >> (hi_bits - 1 - p0)) & 0x1)
                                           : static_cast<uint8_t>((payload_lo >> (WORD_BIT - 1 - (p0 - hi_bits))) & 0x1);
                    const uint8_t b1 = (p1 < hi_bits)
                                           ? static_cast<uint8_t>((payload_hi >> (hi_bits - 1 - p1)) & 0x1)
                                           : static_cast<uint8_t>((payload_lo >> (WORD_BIT - 1 - (p1 - hi_bits))) & 0x1);

                    const uint8_t enc = static_cast<uint8_t>((b0 << 1) | b1);
                    seq.push_back(GlobalVariant::bin2char[enc]);
                }
                Segment segment(node_num, seq, so, sr);
                segment.setSosrStorage(has_sosr, !sosr_uncompressed);
                if (has_sosr)
                    segment.setSnId(sn_id);
                bp_num_ += segment.getLength();
                segments_.push_back(segment);
            }
            else
            {
                throw std::runtime_error("[Load From bGFA File] Invalid segment compact tag in no-id mode.");
            }
        }
        else
        {
            std::vector<uword_t> seq_data(1, 0);

            if ((l_str >> (WORD_BIT - 1)) == 1)
            {
                seq_data[0] = l_str << SEGMENT_COMPACT_PLACEHOLDER;
                l_str = (l_str >> (WORD_BIT - SEGMENT_COMPACT_PLACEHOLDER)) &
                        uword_t((1 << (SEGMENT_COMPACT_PLACEHOLDER - 1)) - 1);
            }
            else
            {
                uint8_t baseinElement = WORD_BIT / 2;
                size_t seq_num =
                    (l_str + baseinElement - 1) / baseinElement;
                seq_data.resize(seq_num);
                if (seq_num > 0)
                {
                    file.read(reinterpret_cast<char *>(seq_data.data()),
                              seq_num * sizeof(uword_t));
                    if (!file)
                    {
                        throw std::runtime_error("[Load From bGFA File] Error reading sequence data in S record.");
                    }
                }
                curr_pos += seq_num;
            }
            Segment segment(node_num, so, sr, l_str, seq_data);
            segment.setSosrStorage(has_sosr, !sosr_uncompressed);
            if (has_sosr)
                segment.setSnId(sn_id);
            bp_num_ += segment.getLength();
            segments_.push_back(segment);
        }

        name_to_id_[std::to_string(node_num)] = seg_num_;
        seg_num_++;
    }

    /**
     * @brief Parse a link line from GFA file
     * @param {string} &line - Line from GFA file
     */
    void parseLinkLine(const std::string &line)
    {
        // Ensure links set sized to current node count
        links_.resizeLinks(seg_num_);
#ifdef INFO_DEBUG
        std::cout << "\tline:" << link_num_ << " Links_size:" << links_.getAllLinks().size() << std::endl;
        if (seg_num_ != links_.getAllLinks().size())
        {
            std::cout << "\tNode_num:" << seg_num_ << " Links_size:" << links_.getAllLinks().size() << std::endl;
            throw std::logic_error("[Parse Link Line] Node number mismatch when parsing L record.");
        }
#endif
        std::istringstream iss(line);
        std::string type, from_name, from_dir, to_name, to_dir, overlap;

        iss >> type >> from_name >> from_dir >> to_name >> to_dir >> overlap;

        uword_t from_id = getSegmentIndex(from_name);
        uword_t to_id = getSegmentIndex(to_name);

        link_num_++;
        links_.addLink(from_id, to_id, from_dir, to_dir);
    }

    /**
     * @brief Parse a binary link record from bGFA file
     * @param {ifstream} &file - Input binary file stream
     * @param {uword_t} &curr_pos - Current position in the file (in uword_t units)
     * @param {uword_t} limit - Limit position for link records
     */
    void parseBinLinkField(std::ifstream &file, uword_t &curr_pos, const uword_t limit)
    {
        uword_t links_num, n_rows;

        file.read(reinterpret_cast<char *>(&links_num), sizeof(links_num));
        file.read(reinterpret_cast<char *>(&n_rows), sizeof(n_rows));

        links_.resizeLinks(n_rows);

        if (seg_num_ != n_rows)
        {
            throw std::logic_error("[Load From bGFA File] Node number mismatch when reading L record.");
        }

        std::vector<uword_t> indptr(n_rows + 1, 0);
        for (int i = 0; i <= n_rows; i++)
        {
            file.read(reinterpret_cast<char *>(&indptr[i]), sizeof(uword_t));

            if (!file)
            {
                throw std::runtime_error("[Load From bGFA File] Error reading L record index data.");
            }
        }

        uword_t curr_link_index = 0;
        for (uword_t i = 0; i < links_num; i++)
        {
            uword_t val;
            file.read(reinterpret_cast<char *>(&val), sizeof(val));
            if (!file)
            {
                throw std::runtime_error("[Load From bGFA File] Error reading L record data.");
            }

            while (i >= indptr[curr_link_index + 1])
            {
                curr_link_index++;
            }
            links_.addLink(curr_link_index, val);
        }

        curr_pos += 2 + links_num + n_rows + 1;
        link_num_ = links_num;
    }

    /**
     * @brief Parse a path line from GFA file
     * @param {string} &line - Line from GFA file
     */
    void parsePathLine(const std::string &line)
    {
        std::istringstream iss(line);
        std::string type, name, path_str, cigar_str;

        iss >> type >> name >> path_str >> cigar_str;
        Path path;
        if (pre_info_ & PRE_INFO_PATHMODE_INC)
            path = Path(path_num_, path_str, name_to_id_, PATHMODE_INCREMENT);
        else
            path = Path(path_num_, path_str, name_to_id_, PATHMODE_DIRECT);
        path_num_++;
        path_size_ += path.getStoreSize();
        paths_.push_back(path);
    }

    /**
     * @brief Parse a binary path record from bGFA file
     * @param {ifstream} &file - Input binary file stream
     * @param {uword_t} &curr_pos - Current position in the file (in uword_t units)
     * @param {uword_t} limit - Limit position for path records
     */
    void parseBinPathField(std::ifstream &file, uword_t &curr_pos, const uword_t limit)
    {
        uword_t n_nodes = 0;
        file.read(reinterpret_cast<char *>(&n_nodes), sizeof(n_nodes));
        if (!file)
        {
            throw std::runtime_error("[Load From bGFA File] Error reading P record header.");
        }
        curr_pos += 1;

        std::vector<uword_t> nodes;
        nodes.reserve(n_nodes);
        for (uword_t i = 0; i < n_nodes; ++i)
        {
            uword_t node_val = 0;
            file.read(reinterpret_cast<char *>(&node_val), sizeof(node_val));
            if (!file)
            {
                throw std::runtime_error("[Load From bGFA File] Error reading P record node value.");
            }
            nodes.emplace_back(node_val);
        }
        curr_pos += n_nodes;

        paths_.emplace_back(Path(nodes));
    }

    /**
     * @brief Parse a walk line from GFA file
     * @param {string} &line - Line from GFA file
     */
    void parseWalkLine(const std::string &line)
    {
        std::istringstream iss(line);
        std::string type, sample_name, seq_name, walk_str;

        uword_t sample_id, hap_id, seq_id;

        iss >> type >> sample_name >> hap_id >> seq_name >> walk_str;

        std::string token;
        while (iss >> token)
        {
            if (!token.empty() && (token[0] == '>' || token[0] == '<'))
            {
                walk_str = token;
            }
        }

        if (w_sample_to_id_.find(sample_name) == w_sample_to_id_.end())
        {
            sample_id = w_sample_to_id_.size();
            w_sample_to_id_[sample_name] = sample_id;
        }
        else
        {
            sample_id = w_sample_to_id_[sample_name];
        }

        if (w_sequence_to_id_.find(seq_name) == w_sequence_to_id_.end())
        {
            seq_id = w_sequence_to_id_.size();
            w_sequence_to_id_[seq_name] = seq_id;
        }
        else
        {
            seq_id = w_sequence_to_id_[seq_name];
        }

        Path path;

        if (pre_info_ & PRE_INFO_PATHMODE_INC)
            path = Path(walk_str, name_to_id_, PATHMODE_INCREMENT);
        else
            path = Path(walk_str, name_to_id_, PATHMODE_DIRECT);
        Walk walk(sample_id, hap_id, seq_id, path);

        walk_num_++;
        walk_size_ += walk.getStoreSize();
        walks_.push_back(walk);
    }

    void parseBinWalkField(std::ifstream &file, uword_t &curr_pos, const uword_t limit)
    {
        uword_t sample_id, hap_id, seq_id, n_nodes;
        file.read(reinterpret_cast<char *>(&sample_id), sizeof(sample_id));
        file.read(reinterpret_cast<char *>(&hap_id), sizeof(hap_id));
        file.read(reinterpret_cast<char *>(&seq_id), sizeof(seq_id));
        file.read(reinterpret_cast<char *>(&n_nodes), sizeof(n_nodes));

        if (!file)
        {
            throw std::runtime_error("[Load From bGFA File] Error reading W record header.");
        }
        curr_pos += 4;

        std::vector<uword_t> nodes;
        nodes.reserve(n_nodes);
        for (uword_t i = 0; i < n_nodes; ++i)
        {
            uword_t node_val = 0;
            file.read(reinterpret_cast<char *>(&node_val), sizeof(node_val));
            if (!file)
            {
                throw std::runtime_error("[Load From bGFA File] Error reading W record node value.");
            }
            nodes.emplace_back(node_val);
        }
        curr_pos += n_nodes;

        Path path(nodes);
        walks_.emplace_back(Walk(sample_id, hap_id, seq_id, path));
    }
};
// Reserve info: last bit indicates section tags presence

bool operator==(const GFA &a, const GFA &b)
{
    // 1) Collect valid segments (length>0) and check count equality
    const auto &segsA = a.getAllSegments();
    const auto &segsB = b.getAllSegments();
    std::vector<uword_t> validA, validB;
    validA.reserve(segsA.size());
    validB.reserve(segsB.size());
    for (uword_t i = 0; i < segsA.size(); ++i)
        if (segsA[i].getLength() > 0)
            validA.push_back(i);
    for (uword_t i = 0; i < segsB.size(); ++i)
        if (segsB[i].getLength() > 0)
            validB.push_back(i);
    if (validA.size() != validB.size())
        return false;

    // 2) Check link count equality
    const auto &linksA = a.getLinksSet();
    const auto &linksB = b.getLinksSet();
    if (linksA.getLinkNum() != linksB.getLinkNum())
        return false;

    // 3) Construct a bijection mapping based on sequences oldA_id -> oldB_id
    std::unordered_map<std::string, uword_t> seq2idB;
    seq2idB.reserve(validB.size());
    for (uword_t idB : validB)
    {
        const std::string s = segsB[idB].getSequenceAsString();
        // Require unique sequences, otherwise stable mapping cannot be established
        if (seq2idB.find(s) != seq2idB.end())
            return false;
        seq2idB.emplace(s, idB);
    }
    std::unordered_map<uword_t, uword_t> mapA2B;
    mapA2B.reserve(validA.size());
    for (uword_t idA : validA)
    {
        const std::string s = segsA[idA].getSequenceAsString();
        auto it = seq2idB.find(s);
        if (it == seq2idB.end())
            return false; // A sequence in A does not exist in B
        mapA2B.emplace(idA, it->second);
    }

    // 4) Compare all edges based on the mapping
    const auto &allA = linksA.getAllLinks();
    const auto &allB = linksB.getAllLinks();

    // Helper: organize B's adjacency sets by from_id for easy lookup
    auto hasEdgeInB = [&](uword_t fromB, uword_t toB, uword_t dirBits) -> bool
    {
        if (fromB >= allB.size())
            return false;
        const auto &setB = allB[fromB];
        uword_t enc = (toB << LINK_DIR_BIT) | (dirBits & 0x03);
        return setB.find(enc) != setB.end();
    };

    for (uword_t fromA = 0; fromA < allA.size(); ++fromA)
    {
        // Only compare outgoing edges of valid nodes; invalid nodes (length==0) are considered non-existent
        auto itFromMap = mapA2B.find(fromA);
        if (itFromMap == mapA2B.end())
        {
            // If this node in A is invalid, it should have no outgoing edges; otherwise not equal
            if (!allA[fromA].empty())
                return false;
            continue;
        }
        uword_t fromB = itFromMap->second;
        const auto &outA = allA[fromA];

        // Each edge in A must have a corresponding edge in B
        for (const auto &valA : outA)
        {
            uword_t toA = (valA >> LINK_DIR_BIT);
            uword_t dirA = (valA & 0x03);
            auto itToMap = mapA2B.find(toA);
            if (itToMap == mapA2B.end())
                return false; // Points to an invalid or unmapped node
            uword_t toB = itToMap->second;
            if (!hasEdgeInB(fromB, toB, dirA))
                return false;
        }

        size_t mappedCountB = 0;
        for (const auto &valB : allB[fromB])
        {
            uword_t toB = (valB >> LINK_DIR_BIT);
            bool existsInA = false;
            for (const auto &kv : mapA2B)
            {
                if (kv.second == toB)
                {
                    // Check if there exists fromA->kv.first (same dir) in A
                    uword_t encA = (kv.first << LINK_DIR_BIT) | (valB & 0x03);
                    // Check if outA contains the corresponding edge
                    if (outA.find(encA) != outA.end())
                    {
                        existsInA = true;
                        break;
                    }
                }
            }
            if (existsInA)
                mappedCountB++;
        }
        if (mappedCountB != outA.size())
            return false;
    }

    return true;
}

#endif // GRAPH_HPP_
