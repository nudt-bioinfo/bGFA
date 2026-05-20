#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>

#include "../bgfa_graph.hpp"
#include "../bgfa_merge.hpp"
#include "../bgfa_alignment.hpp"

static std::vector<std::string> canonicalizeGfa(const std::string &path)
{
    std::ifstream f(path);
    std::vector<std::string> out;
    if (!f.is_open())
        return out;
    // First pass: build name->idx by S order
    std::vector<std::pair<std::string, std::string>> s_lines; // name, seq
    std::vector<std::string> raw_lines;
    raw_lines.reserve(256);
    std::string line;
    while (getline(f, line))
    {
        if (line.empty() || line[0] == '#')
            continue;
        raw_lines.push_back(line);
        if (line[0] == 'S')
        {
            // S <name> <seq>
            std::string type, name, seq;
            std::stringstream ss(line);
            ss >> type >> name >> seq;
            s_lines.emplace_back(name, seq);
        }
    }
    // Build mapping
    std::unordered_map<std::string, int> name2idx;
    name2idx.reserve(s_lines.size() * 2);
    for (size_t i = 0; i < s_lines.size(); ++i)
        name2idx[s_lines[i].first] = static_cast<int>(i);

    // Second pass: emit canonical lines
    for (const auto &ln : raw_lines)
    {
        if (ln[0] == 'S')
        {
            std::string type, name, seq;
            std::stringstream ss(ln);
            ss >> type >> name >> seq;
            int idx = name2idx[name];
            out.push_back("S " + std::to_string(idx) + " " + seq);
        }
        else if (ln[0] == 'L')
        {
            // L <from> <from_dir> <to> <to_dir> <ovl>
            std::string type, from, fdir, to, tdir, ovl;
            std::stringstream ss(ln);
            ss >> type >> from >> fdir >> to >> tdir >> ovl;
            int fi = name2idx[from];
            int ti = name2idx[to];
            out.push_back("L " + std::to_string(fi) + " " + fdir + " " + std::to_string(ti) + " " + tdir + " *");
        }
        else if (ln[0] == 'P')
        {
            // P <id> <nodes> <cigar>
            std::string type, pid, nodes, cigar;
            std::stringstream ss(ln);
            ss >> type >> pid >> nodes >> cigar;
            // rewrite nodes list
            std::stringstream ns(nodes);
            std::string node;
            std::vector<std::string> normNodes;
            while (getline(ns, node, ','))
            {
                if (node.empty())
                    continue;
                char dir = node.back();
                std::string name = node.substr(0, node.size() - 1);
                int ni = name2idx[name];
                normNodes.push_back(std::to_string(ni) + std::string(1, dir));
            }
            std::string joined;
            for (size_t i = 0; i < normNodes.size(); ++i)
            {
                if (i)
                    joined.push_back(',');
                joined += normNodes[i];
            }
            out.push_back("P 0 " + joined + " *");
        }
    }
    return out;
}

static std::string normalize(const std::string &s)
{
    std::string out;
    out.reserve(s.size());
    bool inSpace = false;
    for (char c : s)
    {
        if (c == '\t' || c == ' ')
        {
            if (!inSpace)
            {
                out.push_back(' ');
                inSpace = true;
            }
        }
        else
        {
            out.push_back(c);
            inSpace = false;
        }
    }
    if (!out.empty() && out.front() == ' ')
        out.erase(out.begin());
    if (!out.empty() && out.back() == ' ')
        out.pop_back();
    return out;
}

bool filesEqual(const std::string &a, const std::string &b)
{
    auto ca = canonicalizeGfa(a);
    auto cb = canonicalizeGfa(b);
    if (ca.size() != cb.size())
        return false;
    for (size_t i = 0; i < ca.size(); ++i)
    {
        if (normalize(ca[i]) != normalize(cb[i]))
            return false;
    }
    return true;
}

static std::vector<std::string> readLines(const std::string &path)
{
    std::ifstream f(path);
    std::vector<std::string> out;
    std::string line;
    while (getline(f, line))
    {
        if (line.empty() || line[0] == '#')
            continue;
        out.push_back(line);
    }
    return out;
}

static bool gfaEqualLoose(const std::string &a, const std::string &b)
{
    auto la = readLines(a);
    auto lb = readLines(b);
    if (la.size() != lb.size())
        return false;
    for (size_t i = 0; i < la.size(); ++i)
    {
        if (normalize(la[i]) != normalize(lb[i]))
            return false;
    }
    return true;
}