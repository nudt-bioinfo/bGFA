/*
 * @file /gfatools-bgfa/bgfa.hpp
 * @brief Basic definitions for bGFA project
 * @details This file contains basic type definitions, constants,
 *          and includes necessary for the bGFA project.
 * @author grampart
 * @version V1.0
 * @copyright: Copyright(c) 2025 by grampart, All Rights Reserved.
 * @license: MIT License
 */
#ifndef BGFA_H_
#define BGFA_H_

#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <queue>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <unordered_map>
#include <unordered_set>
#include <unistd.h>
#include <vector>

// Define WORD_BIT if not already defined

#ifndef WORD_BIT
#define WORD_BIT 32 ///< Default to 32 bits if not defined
#endif

#ifndef __WORDSIZE
#define __WORDSIZE WORD_BIT
#endif

#if WORD_BIT == 16
typedef uint16_t uword_t;
#elif WORD_BIT == 32
typedef uint32_t uword_t;
#elif WORD_BIT == 64
typedef uint64_t uword_t;
#endif

#endif // BGFA_H_