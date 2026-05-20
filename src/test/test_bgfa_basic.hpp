#pragma once
#include "test_bgfa.hpp"

static std::string basicDir = "./data/graph/basic";

static std::string makeRepeatACGT(size_t len)
{
    static const std::string motif = "ACGT";
    std::string out;
    out.reserve(len);
    for (size_t i = 0; i < len; ++i)
    {
        out.push_back(motif[i % motif.size()]);
    }
    return out;
}

int test_roundtrip_gfa_bgfa_gfa(const std::string &srcName)
{
    std::string gfaIn = basicDir + "/" + srcName;
    std::string bgfaOut = basicDir + "/" + srcName + ".bgfa";
    std::string gfaOut = basicDir + "/" + srcName + ".roundtrip.gfa";
    GFA gfa(gfaIn);
    gfa.write2Bgfa(bgfaOut);

    GFA gfa2(bgfaOut);
    gfa2.write2GFA(gfaOut);

    if (gfa == gfa2)
    {
        std::cout << "\t[PASS]";
        return 1;
    }
    else
    {
        std::cout << "[Failed] GFA mismatch after roundtrip: " + srcName
                  << std::endl;
        return 0;
    }
}

int test_segments_compact_design()
{
    // sample_segments.gfa contains one short and one long sequence
    std::string src = "g1.gfa";
    std::string bgfaOut = basicDir + "/" + src + ".bgfa";

    GFA gfa(basicDir + "/" + src);
    gfa.write2Bgfa(bgfaOut);

    // Reload and check lengths preserved
    GFA gfa2(bgfaOut);
    const auto &segs = gfa2.getAllSegments();
    if (segs.size() < 2)
    {
        std::cout << "[Failed] Segment count too small " << std::endl;
        return 0;
    }
    bool has4 = false, has32 = false, has3 = false;
    for (const auto &s : segs)
    {
        if (s.getLength() == 4)
            has4 = true;
        if (s.getLength() == 32)
            has32 = true;
        if (s.getLength() == 3)
            has3 = true;
    }
    if (!has4 || !has32)
    {
        std::cout << "[Failed] Missing expected segment lengths (4 and 32) " << std::endl;
        return 0;
    }
    if (!has3)
    {
        std::cout << "[Failed] Missing auxiliary short segment length=3 " << std::endl;
        return 0;
    }

    std::cout << "\t[PASS]";
    return 1;
}

int test_links_design()
{
    std::string src = "g2.gfa";
    GFA gfa(basicDir + "/" + src);
    const auto &links = gfa.getLinksSet();
    if (links.getLinkNum() != 1)
    {
        std::cout << "[Failed] Link count mismatch " << std::endl;
        return 0;
    }
    if (links.getLinkToByIndex(0, 0) != 1)
    {
        std::cout << "[Failed] Target id mismatch " << std::endl;
        return 0;
    }
    if (links.getLinkDirByIndex(0, 0) != 2)
    {
        std::cout << "[Failed] Direction bits mismatch " << std::endl;
        return 0;
    }
    std::cout << "\t[PASS]";
    return 1;
}

int test_paths_design()
{
    std::string src = "g3.gfa";
    GFA gfa(basicDir + "/" + src);
    const auto &paths = gfa.getAllPaths();
    if (paths.size() != 1)
    {
        std::cout << "[Failed] Path count mismatch " << std::endl;
        return 0;
    }
    const auto &nodes = paths[0].getAllNodes();
    if ((nodes[0] >> 1) != 0)
    {
        std::cout << "[Failed] From node id incorrect " << std::endl;
        return 0;
    }
    if ((nodes[1] >> 1) != 2)
    {
        std::cout << "[Failed] End node id incorrect " << std::endl;
        return 0;
    }
    std::cout << "\t[PASS]";
    return 1;
}

int test_path_increment()
{
    // Increment mode nodes layout test

    // Construct a simple NodeMap and LinksSet, then parse in increment mode
    NodeMap id_to_id;
    id_to_id["A"] = 0;
    id_to_id["B"] = 1;
    id_to_id["C"] = 2;
    id_to_id["D"] = 3;
    id_to_id["E"] = 4;
    id_to_id["F"] = 5;

    std::string path_str = "A+,B+,C+,D+,E+,F-";
    std::vector<uword_t> res = {0x2, 0x14, 0x7, 0xbf, 0x1010101};
    Path p(0, path_str, id_to_id, PATHMODE_INCREMENT);
    const auto &nodes = p.getAllNodes();
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        if (nodes[i] != res[i])
        {
            std::cout << "[Failed] Increment mode node encoding mismatch at index " << i << std::endl;
            return 0;
        }
    }

    std::cout << p.toGfaTxt();

    return 1;
}

int test_segment_no_id_encoding_00_10_11()
{
    const std::string gfaIn = "./temp/test_noid_001011.gfa";
    const std::string bgfaOut = "./temp/test_noid_001011.bgfa";

    const std::string seq10 = makeRepeatACGT(5);  // <=13, tag 10
    const std::string seq11 = makeRepeatACGT(20); // 14..28, tag 11
    const std::string seq00 = makeRepeatACGT(40); // >28, tag 00

    {
        std::ofstream ofs(gfaIn);
        if (!ofs.is_open())
        {
            std::cout << "[Failed] Cannot create test gfa file " << gfaIn << std::endl;
            return 0;
        }
        ofs << "H\tVN:Z:1.0\n";
        ofs << "S\t1\t" << seq10 << "\n";
        ofs << "S\t2\t" << seq11 << "\n";
        ofs << "S\t3\t" << seq00 << "\n";
    }

    GFA gfa(gfaIn);
    gfa.write2Bgfa(bgfaOut, false);

    {
        std::ifstream ifs(bgfaOut, std::ios::binary);
        if (!ifs.is_open())
        {
            std::cout << "[Failed] Cannot open bgfa output " << bgfaOut << std::endl;
            return 0;
        }

        uword_t pre_info = 0, seg_size = 0, link_size = 0, path_size = 0, walk_size = 0;
        ifs.read(reinterpret_cast<char *>(&pre_info), sizeof(uword_t));
        ifs.read(reinterpret_cast<char *>(&seg_size), sizeof(uword_t));
        ifs.read(reinterpret_cast<char *>(&link_size), sizeof(uword_t));
        ifs.read(reinterpret_cast<char *>(&path_size), sizeof(uword_t));
        ifs.read(reinterpret_cast<char *>(&walk_size), sizeof(uword_t));
        if (!ifs)
        {
            std::cout << "[Failed] Cannot read bgfa header" << std::endl;
            return 0;
        }

        if ((pre_info & PRE_INFO_SEGMENT_NO_ID) == 0)
        {
            std::cout << "[Failed] PRE_INFO_SEGMENT_NO_ID not set" << std::endl;
            return 0;
        }

        uword_t sosr = 0, l1 = 0, l2 = 0, l3 = 0, payload_lo = 0, seqw = 0;
        ifs.read(reinterpret_cast<char *>(&sosr), sizeof(uword_t));
        ifs.read(reinterpret_cast<char *>(&l1), sizeof(uword_t));
        if ((l1 >> (WORD_BIT - 2)) != 2)
        {
            std::cout << "[Failed] Expected tag 10 for first segment" << std::endl;
            return 0;
        }

        ifs.read(reinterpret_cast<char *>(&sosr), sizeof(uword_t));
        ifs.read(reinterpret_cast<char *>(&l2), sizeof(uword_t));
        if ((l2 >> (WORD_BIT - 2)) != 3)
        {
            std::cout << "[Failed] Expected tag 11 for second segment" << std::endl;
            return 0;
        }
        ifs.read(reinterpret_cast<char *>(&payload_lo), sizeof(uword_t));

        ifs.read(reinterpret_cast<char *>(&sosr), sizeof(uword_t));
        ifs.read(reinterpret_cast<char *>(&l3), sizeof(uword_t));
        if (l3 != 40)
        {
            std::cout << "[Failed] Expected plain length(00 branch)=40, got " << l3 << std::endl;
            return 0;
        }
        for (int i = 0; i < 3; ++i)
            ifs.read(reinterpret_cast<char *>(&seqw), sizeof(uword_t));

        if (!ifs)
        {
            std::cout << "[Failed] Binary payload parsing failed" << std::endl;
            return 0;
        }
    }

    GFA gfa2(bgfaOut);
    const auto &segs = gfa2.getAllSegments();
    if (segs.size() != 3)
    {
        std::cout << "[Failed] Segment count mismatch after reload" << std::endl;
        return 0;
    }

    if (segs[0].getId() != 0 || segs[1].getId() != 1 || segs[2].getId() != 2)
    {
        std::cout << "[Failed] Auto-assigned segment ids are not sequential" << std::endl;
        return 0;
    }

    if (segs[0].getSequenceAsString() != seq10 ||
        segs[1].getSequenceAsString() != seq11 ||
        segs[2].getSequenceAsString() != seq00)
    {
        std::cout << "[Failed] Sequence mismatch after no-id roundtrip" << std::endl;
        return 0;
    }

    std::cout << "\t[PASS]";
    return 1;
}

void test_bgfa_basic(int &pass_num, int &test_num)
{
    test_num += 8;

    std::cout << "Running BGFA Basic Graph Tests..." << std::endl;

    pass_num += test_roundtrip_gfa_bgfa_gfa("g1.gfa");
    std::cout << " Basic BGFA Graph Tests[1]" << std::endl;

    pass_num += test_roundtrip_gfa_bgfa_gfa("g2.gfa");
    std::cout << " Basic BGFA Graph Tests[2]" << std::endl;

    pass_num += test_roundtrip_gfa_bgfa_gfa("g3.gfa");
    std::cout << " Basic BGFA Graph Tests[3]" << std::endl;

    pass_num += test_segments_compact_design();
    std::cout << " Segment Compact Test" << std::endl;

    pass_num += test_links_design();
    std::cout << " Basic Link Test" << std::endl;

    pass_num += test_paths_design();
    std::cout << " Basic Path Test" << std::endl;

    pass_num += test_path_increment();
    std::cout << " Path Increment Mode Test" << std::endl;

    pass_num += test_segment_no_id_encoding_00_10_11();
    std::cout << " Segment No-ID Encoding 00/10/11 Test" << std::endl;

    return;
}
