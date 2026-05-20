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

/// Placeholder size, used for representing link direction
#define LINK_DIR_BIT 2
/// Placeholder size, used for flag/length in compact segme nt representation
#define SEGMENT_COMPACT_PLACEHOLDER 6
/// Placeholder size, used for offset/rank in compact segment representation
#define SOSR_OFFSET_SHIFT 6
/// Mask for rank in sosr
#define SR_MASK 0x3F

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

/**
 * @brief Segment class representing a graph segment
 * @details The equality comparison for segments does not follow the same rules
 *          as the greater-than and less-than comparisons. Two segments are
 *          considered equal if their sequences are identical, regardless of
 *          their IDs or other attributes. A segment is considered less than
 *          another if its dis is smaller, and then by id.
 */
class BGFA_Segment
{
private:
    uword_t id_;                    ///< Segment ID
    uword_t dis_;                   ///< Segment distance, not stored
    uword_t sosr_;                  ///< Segment offset and rank (low 6 bits for rank, high bits for offset)
    uword_t length_;                ///< Segment length(SL:i), not stored
    std::vector<uword_t> sequence_; ///< Binary encoded sequence

    uword_t store_size_; ///< Size in storage (in uword_t units)

public:
    /**
     * @brief Default constructor initializes an empty segment
     */
    BGFA_Segment() : length_(0), store_size_(0) {}

    /**
     * @brief Constructor with id, dis, length, and binary sequence
     * @param {uword_t} id - Segment ID
     * @param {uword_t} sosr - Segment offset and rank
     * @param {uword_t} length - Segment length
     * @param {vector<uword_t>} &sequence - Binary encoded sequence
     */
    explicit BGFA_Segment(const uword_t id, uword_t sosr, uword_t length,
                          const std::vector<uword_t> &sequence)
        : id_(id), length_(length), sequence_(sequence), sosr_(sosr)
    {
        setStoreSize();
    }

    /**
     * @brief Constructor with id, length, and binary sequence (default: 0)
     * @param {uword_t} id - Segment ID
     * @param {uword_t} length - Segment length
     * @param {vector<uword_t>} &sequence - Binary encoded sequence
     */
    explicit BGFA_Segment(const uword_t id, uword_t length,
                          const std::vector<uword_t> &sequence)
        : id_(id), length_(length), sequence_(sequence), dis_(0), sosr_(0)
    {
        setStoreSize();
    }

    /**
     * @brief Constructor with id and string sequence (dis default: 0)
     * @param {uword_t} id - Segment ID
     * @param {string} &sequence_str - Segment sequence string
     */
    explicit BGFA_Segment(const uword_t id, const std::string &sequence_str)
        : id_(id), length_(sequence_str.length()), dis_(0), sosr_(0)
    {
        stringToBinary(sequence_str);
        setStoreSize();
    }

    explicit BGFA_Segment(const uword_t id, const std::string &sequence_str, const uword_t offset, const uword_t rank)
        : id_(id), length_(sequence_str.length()), dis_(0)
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
    uword_t getSosr() const { return sosr_; }

    /**
     * @brief Get the Offset of the Segment
     * @return {uword_t} Segment offset
     */
    uword_t getOffset() const { return sosr_ >> SOSR_OFFSET_SHIFT; }

    /**
     * @brief Get the Rank of the Segment
     * @return {uword_t} Segment rank
     */
    uword_t getRank() const { return sosr_ & SR_MASK; }

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
        sosr_ = (offset << SOSR_OFFSET_SHIFT) | (rank & SR_MASK);
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
        oss << "\tSZ:Z:" << getRank();

        return oss.str();
    }

    /**
     * @brief Convert the Segment to bGFA binary format
     * @param {ofstream} &output_file - Output binary file stream
     * @param {uword_t} &s_offset - Current offset in the binary file
     */
    void toBgfaBin(std::ofstream &output_file) const
    {
        if (!output_file.is_open())
        {
            throw std::runtime_error(
                "[Write to bGFA File] Output file is not open");
        }

        uint8_t type_flag = 1; // Type flag for segment

        // output_file.write(reinterpret_cast<const char *>(&type_flag), sizeof(uword_t));
        output_file.write(reinterpret_cast<const char *>(&id_),
                          sizeof(uword_t));
        output_file.write(reinterpret_cast<const char *>(&sosr_),
                          sizeof(uword_t));

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
            uword_t val = (((1 << 5) | length_) << (WORD_BIT - 6)) |
                          (sequence_[0] >> 6);
            output_file.write(reinterpret_cast<const char *>(&val),
                              sizeof(uword_t));
        }
    }

    /**
     * @brief Equality operator for Segment
     * @details Two segments are equal if their sequences are same.
     * @param {Segment} &other - Other segment to compare
     * @return {bool} True if segments are equal, false otherwise
     */
    bool operator==(const BGFA_Segment &other) const
    {
        return sequence_ == other.sequence_;
    }

    /**
     * @brief Inequality operator for Segment
     * @param {Segment} &other - Other segment to compare
     * @return {bool} True if segments are not equal, false otherwise
     */
    bool operator!=(const BGFA_Segment &other) const
    {
        return !(*this == other);
    }

    /**
     * @brief Less-than operator for Segment (based on distance, then ID)
     * @param {Segment} &other - Other segment to compare
     * @return {bool}
     */
    bool operator<(const BGFA_Segment &other) const
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
    bool operator>(const BGFA_Segment &other) const { return other < *this; }

    /**
     * @brief Less-than-or-equal-to operator for Segment
     * @param {Segment} &other - Other segment to compare
     * @return {bool}
     */
    bool operator<=(const BGFA_Segment &other) const { return !(other < *this); }

    /**
     * @brief Greater-than-or-equal-to operator for Segment
     * @param {Segment} &other - Other segment to compare
     * @return {bool}
     */
    bool operator>=(const BGFA_Segment &other) const { return !(*this < other); }

    /**
     * @brief Output stream operator for Segment
     * @param {ostream} &os - Output stream
     * @param {Segment} &segment - Segment to output
     * @return {ostream&} Reference to the output stream
     */
    friend std::ostream &operator<<(std::ostream &os, const BGFA_Segment &segment)
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

        os << "\n\tSO: " << (segment.sosr_ >> SOSR_OFFSET_SHIFT)
           << "\tSR: " << (segment.sosr_ & SR_MASK);

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

        store_size_ += 1; // for sosr

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
class BGFA_LinksSet
{
private:
    uword_t link_num_;                               ///< Number of links
    uword_t store_size_;                             ///< Size in storage (in uword_t units)
    std::vector<std::unordered_set<uword_t>> links_; ///< Adjacency set representation of links, links[from_id] contains encoded (to_id<<2)|dir

public:
    /**
     * @brief Default constructor initializes an empty LinksSet
     */
    BGFA_LinksSet() : link_num_(0), store_size_(0) {}

    /**
     * @brief Constructor that builds LinksSet from given data and indptr vectors
     * @param {uword_t} links_num - Number of links
     * @param {uword_t} nodes_num - Number of nodes
     * @param {vector<uword_t>} &data - Data vector containing link information
     * @param {vector<uword_t>} &indptr - Indptr vector for indexing
     */
    explicit BGFA_LinksSet(const uword_t links_num, const uword_t nodes_num,
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
    friend std::ostream &operator<<(std::ostream &os, const BGFA_LinksSet &links_set)
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

class BGFA_Path
{
private:
    uword_t id_;
    std::vector<uword_t> nodes_; // (id<<1) | ('+' ? 1 : 0)

    uword_t store_size_; ///< Size in storage (in uword_t units)

public:
    BGFA_Path() = default;
    explicit BGFA_Path(const std::vector<uword_t> &nodes) : nodes_(nodes)
    {
        setStoreSize();
    }

    explicit BGFA_Path(uword_t id, const std::vector<uword_t> &nodes)
        : id_(id), nodes_(nodes)
    {
        setStoreSize();
    }

    explicit BGFA_Path(uword_t path_id, std::string path_str, BGFA_LinksSet &links_set, NodeMap &id_to_id) : id_(path_id)
    {
        std::istringstream iss(path_str);
        std::string node_str;
        uword_t node, pre_node_id, node_id;
        char pre_dir, dir;

        // from id
        std::getline(iss, node_str, ',');
        if (node_str.empty())
        {
            throw std::logic_error(
                "[Create Paths Info] Invalid path string: " + path_str);
        }

        pre_dir = node_str.back();
        pre_node_id = id_to_id[node_str.substr(0, node_str.length() - 1)];
        node = (pre_node_id << 1) | (GlobalVariant::dir2bin[pre_dir]);
        nodes_.push_back(node);
        nodes_.push_back(node);

        std::getline(iss, node_str, ',');
        while (1)
        {
            // if (node_str.empty())

            dir = node_str.back();
            node_id = id_to_id[node_str.substr(0, node_str.length() - 1)];
            node = (node_id << 1) | (GlobalVariant::dir2bin[dir]);

            // if (isNotUniqueSuccessor(links_set, pre_node_id,
            //                          node_id, pre_dir, dir))
            // {
            //     nodes_.push_back(node);
            // }

            nodes_.push_back(node);

            pre_node_id = node_id;
            pre_dir = dir;
            if (!std::getline(iss, node_str, ','))
            {
                nodes_[1] = node;
                break;
            }
        }

        setStoreSize();
    }

    explicit BGFA_Path(std::string walk_str, BGFA_LinksSet &links_set, NodeMap &id_to_id) : id_(0)
    {
        // Parse walk_str format: (>|<)id(>|<)id...(>|<)id
        auto is_space = [](char c) -> bool
        { return c == ' ' || c == '\t' || c == '\r' || c == '\n'; };

        size_t pos = 0;
        auto read_next = [&]() -> bool
        {
            // skip whitespaces
            while (pos < walk_str.size() && is_space(walk_str[pos]))
                ++pos;
            if (pos >= walk_str.size())
                return false;

            char dir = walk_str[pos];
            if (dir != '>' && dir != '<')
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
            uword_t dir_bit = (dir == '>') ? 1 : 0; // '>' corresponds to '+', '<' to '-'
            uword_t node = (node_id << 1) | dir_bit;

            if (nodes_.empty())
            {
                nodes_.push_back(node); // from
                nodes_.push_back(node); // placeholder end
            }
            else
            {
                nodes_.push_back(node);
            }
            return true;
        };

        while (read_next())
            ;

        if (!nodes_.empty())
            nodes_[1] = nodes_.back();

        setStoreSize();
    }

    // uword_t getLength() const { return nodes_.size(); }

    /**
     * @brief Get the from node of the BGFA_Path
     * @return {uword_t} From node
     */
    uword_t getFromNode() const
    {
        return nodes_[0];
    }

    /**
     * @brief Get the end node of the BGFA_Path
     * @return {uword_t} End node
     */
    uword_t getEndNode() const { return nodes_[1]; }

    /**
     * @brief Get all nodes in the BGFA_Path
     * @return {const vector<uword_t>&} Vector of all nodes
     */
    const std::vector<uword_t> &getAllNodes() const { return nodes_; }

    /**
     * @brief Get the storage size of the BGFA_Path in uword_t units
     * @return {uword_t} Storage size
     */
    uword_t getStoreSize()
    {
        setStoreSize();
        return store_size_;
    }

    /**
     * @brief Convert the BGFA_Path to GFA text format
     * @param {LinksSet} &links_set - LinksSet for reference
     * @return {string} GFA text representation of the path
     */
    std::string toGfaTxt(const BGFA_LinksSet &links_set) const
    {
        std::ostringstream oss;

        // uword_t path_id =

        oss << "P\t" << 0 << "\t" << (nodes_[0] >> 1)
            << GlobalVariant::bin2dir[nodes_[0] & 0x01] << ",";

        uword_t prev_id = nodes_[0] >> 1;
        uword_t path_index = 2;

        while (prev_id != (nodes_[1] >> 1))
        {
            if (links_set.getLinksFromNode(prev_id).empty())
                throw std::logic_error(
                    "[BGFA_Path to GFA Text] No outgoing links from node: " +
                    std::to_string(prev_id));
            else if (links_set.getLinksFromNode(prev_id).size() == 1)
            {
                uword_t next_id = links_set.getLinkToByIndex(prev_id, 0);
                char next_dir = GlobalVariant::bin2dir[links_set.getLinkDirByIndex(prev_id, 0) & 0x01];

                oss << next_id << next_dir;
                prev_id = next_id;
            }
            else
            {
                oss << (nodes_[path_index] >> 1)
                    << GlobalVariant::bin2dir[nodes_[path_index] & 0x01];
                prev_id = nodes_[path_index] >> 1;
                path_index++;
            }

            if (prev_id != (nodes_[1] >> 1))
                oss << ",";
        }
        oss << " *";
        return oss.str();
    }

    /**
     * @brief Convert the BGFA_Path to bGFA binary format
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
     * @brief Output stream operator for BGFA_Path
     * @param {ostream} &os - Output stream
     * @param {BGFA_Path} &path - BGFA_Path to output
     * @return {ostream&} Reference to the output stream
     */
    friend std::ostream &operator<<(std::ostream &os, const BGFA_Path &path)
    {
        os << "P";

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
    bool isNotUniqueSuccessor(BGFA_LinksSet &links_set,
                              uword_t from_id, uword_t to_id,
                              char from_dir, char to_dir) const
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

    /**
     * @brief Set the storage size of the BGFA_Path in uword_t units
     */
    void setStoreSize()
    {
        store_size_ = 1;              // for n_nodes
        store_size_ += nodes_.size(); // for nodes
    }
};

class BGFA_Walk
{
private:
    BGFA_Path path_;
    uword_t sample_id_;    ///< Sample ID
    uword_t haplotype_id_; ///< Haplotype ID
    uword_t sequence_id_;  ///< Sequence ID
    uword_t store_size_;   ///< Size in storage (in uword_t units)

public:
    BGFA_Walk() : sample_id_(0), haplotype_id_(0), sequence_id_(0), store_size_(0)
    {
        setStoreSize();
    }

    explicit BGFA_Walk(const BGFA_Path &path)
        : path_(path), sample_id_(0), haplotype_id_(0), sequence_id_(0), store_size_(0)
    {
        setStoreSize();
    }

    explicit BGFA_Walk(uword_t sample_id, uword_t haplotype_id, uword_t sequence_id, const BGFA_Path &path)
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

    /**
     * @brief Get the underlying BGFA_Path of the Walk
     * @return {const BGFA_Path&} BGFA_Path object
     */
    const BGFA_Path &getPath() const { return path_; }

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
    std::string toGfaTxt(const BGFA_LinksSet & /*links_set*/) const
    {
        std::ostringstream oss;
        const auto &nodes = path_.getAllNodes();
        if (nodes.size() < 2)
        {
            throw std::logic_error("[Walk to GFA Text] Invalid path nodes in walk");
        }

        // Header: sample and haplotype
        oss << "W\t" << sample_id_ << "\t" << haplotype_id_ << "\t";

        // Print nodes with direction symbol before id, separated by spaces
        auto dir_symbol = [](uword_t node) -> char
        { return (node & 0x01) ? '>' : '<'; };

        // start node
        oss << dir_symbol(nodes[0]) << (nodes[0] >> 1);

        // intermediate choices (skip nodes[1] which stores end placeholder)
        for (size_t i = 2; i < nodes.size(); ++i)
        {
            oss << dir_symbol(nodes[i]) << (nodes[i] >> 1);
        }

        // end marker (kept consistent with BGFA_Path style)
        // oss << " *";

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
    friend std::ostream &operator<<(std::ostream &os, const BGFA_Walk &walk)
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
class BGFA_GFA
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

    std::vector<BGFA_Segment> segments_; ///< List of segments
    BGFA_LinksSet links_;                ///< Set of links
    std::vector<BGFA_Path> paths_;       ///< List of paths
    std::vector<BGFA_Walk> walks_;       ///< List of walks
    NodeMap name_to_id_;                 ///< Mapping from segment names to IDs
    NodeMap w_sample_to_id_;             ///< Mapping from sample names to IDs
    NodeMap w_sequence_to_id_;           ///< Mapping from sample names to IDs

public:
    /**
     * @brief Construct a new GFA object
     */
    BGFA_GFA() : seg_num_(0), bp_num_(0), link_num_(0), path_num_(0), walk_num_(0),
                 seg_size_(0), link_size_(0), path_size_(0), walk_size_(0)
    {
    }

    /**
     * @brief Loads GFA file (supports .gfa and .bgfa formats)
     * @param {string} &filename - BGFA_Path to the GFA file
     */
    explicit BGFA_GFA(const std::string &filename) : BGFA_GFA()
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
    const BGFA_Segment &getSegmentByIndex(uword_t index) const
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
    const std::vector<BGFA_Segment> &getAllSegments() const { return segments_; }

    /**
     * @brief Find segments by distance within a specified percentage difference
     * @param {uword_t} dis - Target distance
     * @param {uword_t} diff_precent - Allowed percentage difference
     * @return {vector<Segment>} Vector of matching segments (returned by value)
     */
    std::vector<BGFA_Segment> findSegmentsByDistance(uword_t dis, uword_t diff_precent)
    {
        std::vector<BGFA_Segment> result_segments;
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
    const BGFA_LinksSet &getLinksSet() const { return links_; }

    /**
     * @brief Get all paths in the GFA
     * @details It returns a constant reference to the vector of BGFA_Path objects
     * @return {const vector<BGFA_Path>&} Vector of all paths
     */
    const std::vector<BGFA_Path> &getAllPaths() const { return paths_; }

    /**
     * @brief Get all walks in the GFA
     * @details It returns a constant reference to the vector of Walk objects
     * @return {const vector<Walk>&} Vector of all walks
     */
    const std::vector<BGFA_Walk> &getAllWalks() const { return walks_; }

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
    uword_t addSegment(BGFA_Segment segment)
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
     * @param {BGFA_Path} path - BGFA_Path object to add
     */
    void addPath(BGFA_Path path)
    {
        paths_.push_back(path);

        path_num_++;
        path_size_ += path.getStoreSize();
    }

    /**
     * @brief Add a new walk to the GFA
     * @param {Walk} walk - Walk object to add
     */
    void addWalk(BGFA_Walk walk)
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
        std::vector<BGFA_Segment> new_segments;
        new_segments.reserve(segments_.size());
        name_to_id_.clear();

        for (uword_t old_id = 0; old_id < segments_.size(); ++old_id)
        {
            const BGFA_Segment &s = segments_[old_id];
            if (s.getLength() == 0)
                continue; // skip deleted
            BGFA_Segment copy(0, s.getSequenceAsString());
            uword_t new_id = new_segments.size();
            copy.setId(new_id);
            copy.setDis(s.getDis());
            new_segments.push_back(copy);
            id_map[old_id] = new_id;
            name_to_id_[std::to_string(new_id)] = new_id;
        }

        // Rebuild links with new IDs
        BGFA_LinksSet new_links;
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
    uword_t isSegmentInGraph(BGFA_Segment seg)
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
    friend std::ostream &operator<<(std::ostream &os, const BGFA_GFA &gfa)
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
            file << path.toGfaTxt(links_) << std::endl;
        }

        for (const auto &walk : walks_)
        {
            file << walk.toGfaTxt(links_) << std::endl;
        }

        file.close();
    }

    /**
     * @brief Write the GFA graph to a bGFA binary format file
     * @param {string} &filename - Output filename
     */
    void write2Bgfa(const std::string &filename)
    {
        std::ofstream ofile(filename, std::ios::binary);
        if (!ofile.is_open())
        {
            throw std::runtime_error(
                "[Write to bGFA File] Failed to open file for writing: " +
                filename);
        }

        setPreInfo(ofile);

        seg_size_ = 0;
        link_size_ = 0;
        path_size_ = 0;
        walk_size_ = 0;

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
            segment.toBgfaBin(ofile);
            seg_size_ += segment.getStoreSize();
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
        file.flush();
        file.close();
    }

private:
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

        link_size_ = links_.getStoreSize();

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
        std::ifstream file(filename);
        uword_t pre_info, s_offset, l_offset, p_offset;
        uword_t curr_pos = 0;

        if (!file.is_open())
        {
            throw std::runtime_error(
                "[Load From bGFA File] Failed to open file: " + filename);
        }

        file.read(reinterpret_cast<char *>(&pre_info), sizeof(pre_info));
        file.read(reinterpret_cast<char *>(&seg_size_), sizeof(seg_size_));
        file.read(reinterpret_cast<char *>(&link_size_), sizeof(link_size_));
        file.read(reinterpret_cast<char *>(&path_size_), sizeof(path_size_));
        file.read(reinterpret_cast<char *>(&walk_size_), sizeof(walk_size_));

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
            // BGFA_Path record
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
        uword_t so = 0, sr = 0;

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
        }

        name_to_id_[id] = seg_num_;
        BGFA_Segment segment(seg_num_, sequence, so, sr);

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
        uword_t node_num, sosr, l_str;
        file.read(reinterpret_cast<char *>(&node_num), sizeof(node_num));
        file.read(reinterpret_cast<char *>(&sosr), sizeof(sosr));
        file.read(reinterpret_cast<char *>(&l_str), sizeof(l_str));

        if (!file)
        {
            throw std::runtime_error(
                "[Load From bGFA File] Error reading S record header.");
        }
        curr_pos += 3;

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
        segments_.push_back(BGFA_Segment(node_num, sosr, l_str, seq_data));
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
        links_.resizeLinks(seg_num_);

        uword_t links_num, val_num;

        file.read(reinterpret_cast<char *>(&links_num), sizeof(links_num));
        file.read(reinterpret_cast<char *>(&val_num), sizeof(val_num));

        if (seg_num_ != val_num)
        {
            throw std::logic_error("[Load From bGFA File] Node number mismatch when reading L record.");
        }

        std::vector<uword_t> indptr(val_num + 1, 0);
        for (int i = 0; i <= val_num; i++)
        {
            file.read(reinterpret_cast<char *>(&indptr[i]), sizeof(uword_t));

            // std::cout << indptr[i] << " ";

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

        curr_pos += 2 + links_num + val_num + 1;
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

        BGFA_Path path(path_num_, path_str, links_, name_to_id_);

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

        paths_.emplace_back(BGFA_Path(nodes));
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

        BGFA_Path path(walk_str, links_, name_to_id_);
        BGFA_Walk walk(sample_id, hap_id, seq_id, path);

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

        BGFA_Path path(nodes);
        walks_.emplace_back(BGFA_Walk(sample_id, hap_id, seq_id, path));
    }

    /**
     * @brief Set pre-info in the bGFA binary file
     * @param {ofstream} &output_file - Output binary file stream
     */
    void setPreInfo(std::ofstream &ofile) const
    {
        // TODO: xxx
        uword_t pre_info = 0;
        ofile.write(reinterpret_cast<const char *>(&pre_info),
                    sizeof(uword_t));
    }
};

bool operator==(const BGFA_GFA &a, const BGFA_GFA &b)
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
