/*
 * @Author: grampart guo018248@163.com
 * @Date: 2025-02-26 15:54:49
 * @LastEditors: grampart guo018248@163.com
 * @LastEditTime: 2025-11-12 16:53:44
 * @FilePath: /bGFA/src/compress.hpp
 * @Description: It Is abandoned
 */
#ifndef COMPRESS_HPP_
#define COMPRESS_HPP_

#include "bgfa.hpp"
#include "bgfa_graph.hpp"
#include <zlib.h>

void zlib_compress_binfile(const std::string &input_path, const std::string &output_path)
{
    // 读取二进制文件
    std::ifstream input_file(input_path, std::ios::binary | std::ios::ate);
    if (!input_file)
    {
        std::cerr << "Failed to open input file: " << input_path << std::endl;
        return;
    }

    std::streamsize file_size = input_file.tellg();
    input_file.seekg(0, std::ios::beg);

    std::vector<char> buffer(file_size);
    if (!input_file.read(buffer.data(), file_size))
    {
        std::cerr << "Failed to read input file: " << input_path << std::endl;
        return;
    }
    input_file.close();

    // 计算压缩后数据所需的大小
    uLongf compressed_size = compressBound(file_size);
    std::vector<Bytef> compressed_data(compressed_size);

    // 进行压缩
    int result = compress(compressed_data.data(), &compressed_size, reinterpret_cast<const Bytef *>(buffer.data()), file_size);
    if (result != Z_OK)
    {
        std::cerr << "Failed to compress data, error code: " << result << std::endl;
        return;
    }

    // 将压缩数据写入文件
    std::ofstream output_file(output_path, std::ios::binary);
    if (!output_file)
    {
        std::cerr << "Failed to open output file: " << output_path << std::endl;
        return;
    }

    output_file.write(reinterpret_cast<const char *>(compressed_data.data()), compressed_size);
    output_file.close();

    std::cout << "Compression successful! Compressed file: " << output_path << std::endl;
    return;
}
#endif