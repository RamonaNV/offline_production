/*The MIT License (MIT)

Copyright (c) 2020, Hendrik Schwanekaml, hschwanekamp@nvidia.com

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGSEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

/* 
    contains settings for the new optimizations
*/

#ifndef SETTINGS_CUH
#define SETTINGS_CUH

// cuda launch params
constexpr int NTHREADS_PER_BLOCK = 512;
constexpr int NBLOCKS_PER_SM = 4;

// wlen lut
constexpr int WLEN_LUT_SIZE = 1024;

// z offset lut
constexpr inline float calcNR(float x, float y) {return -7.0710678119e-01f * x + -7.0710678119e-01f * y;}

constexpr int ZOLUT_NUM_ENTRIES_NR = 128;
constexpr float ZOLUT_MIN_ENTRY_NR = calcNR(800,800);
constexpr float ZOLUT_MAX_ENTRY_NR = calcNR(-800,-800);
constexpr float ZOLUT_SPACING_NR =  (ZOLUT_MAX_ENTRY_NR - ZOLUT_MIN_ENTRY_NR) / float(ZOLUT_NUM_ENTRIES_NR);

constexpr int ZOLUT_NUM_ENTRIES_Z = 256; 
constexpr float ZOLUT_MIN_ENTRY_Z = -800;
constexpr float ZOLUT_MAX_ENTRY_Z = 800;
constexpr float ZOLUT_SPACING_Z = (ZOLUT_MAX_ENTRY_Z - ZOLUT_MIN_ENTRY_Z) / float(ZOLUT_NUM_ENTRIES_Z);

constexpr int ZOLUT_SIZE = ZOLUT_NUM_ENTRIES_NR * ZOLUT_NUM_ENTRIES_Z;

// shared memory usage for look up tables
#define SHARED_WLEN
#define SHARED_ICE_PROPERTIES 
#define SHARED_WLEN_BIAS
#define SHARED_NUM_INDEX_STRING_SET
#define SHARED_COLLISION_GRID_DATA
#define SHARED_STRING_DATA
#define SHARED_STRING_POSITIONS

// blocked random numbers
// #define BLOCK_RANDOM_NUMBERS_SCATTERING
// #define BLOCK_RANDOM_NUMBERS_PROPAGATION

#endif