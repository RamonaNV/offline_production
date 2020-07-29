/*The MIT License (MIT)

Copyright (c) 2020, Hendrik Schwanekamp, hschwanekamp@nvidia.com

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
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef regressionTesting_H
#define regressionTesting_H

#define REPRODUCEABLE_RNG // komment that in / out to enable the same random number to be used for every step

#define REP_RNG_SETS 16384
#define REP_RNG_NUMS_PER_SET 512

#include <vector>
#include <random>

//#define REPRODUCEABLE_RNG // komment that in / out to enable the same random number to be used for every step

// instructions:
// if REPRODUCEABLE_RNG is defined, code should:
// generate numbers using genReproduceableRandomNumbers() below and upload to the device
// have variables "set" and "numInSet" 
// every photon should use a different set calcuated as "set = (stepId+photonId)%REP_RNG_SETS" and start with "numInSet =0" 
// every time a rng number is generated "numInSet = (numInSet+1)%REP_RNG_NUMS_PER_SET"
// get the random number from the array returned by genReproduceableRandomNumbers() at position (set * REP_RNG_NUMS_PER_SET + numInSet)
// that should result in compareable results across multiple runs of all variations of the code

//#define REP_RNG_SETS 16384
//#define REP_RNG_NUMS_PER_SET 512
inline std::vector<float> genReproduceableRandomNumbers();

// generate some random numbers, will generate the same number every time the function is called 
inline std::vector<float> genReproduceableRandomNumbers()
{
    const int numbersTotal = REP_RNG_SETS * REP_RNG_NUMS_PER_SET;
    const int seed = 23072020;
    std::vector<float> reproduceableRandomNumbers(numbersTotal);
    std::default_random_engine rng(seed);
    std::uniform_real_distribution<float> dist(0,1);
    for(int i=0; i<numbersTotal;++i) {
        reproduceableRandomNumbers[i] = dist(rng);
    }
    return reproduceableRandomNumbers;
}



  #endif