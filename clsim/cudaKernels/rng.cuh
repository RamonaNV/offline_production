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
    device random number generation and generator setup on host
    using an multiply with carry rng 
*/

#ifndef RNG_CUH
#define RNG_CUH

// includes
// ------------------
#include <cuda.h>
#include <cuda_runtime_api.h>
// ------------------

/**
 * @brief random number generation on the device using the multiply with carry rng.
 *      Class should live in local memory.
 *      Initialize with a prime number and a state (from global memory), 
 *      then use the different functions to generate random numbers.
 *      Use getState() to store the current state for later use.
 */
class MultipleWithCarryRng 
{
public:
    __device__ MultipleWithCarryRng(uint64_t x, uint32_t a) : m_x(x), m_a(a) {}
    __device__ uint64_t getState() {return m_x;} // returns state (x) to store it back into global memory

    // generates 4 byte of random data
    __device__ __forceinline__ uint32_t rand() 
    { 
        m_x = (m_x & 0xffffffffull) * m_a + (m_x >> 32);
        return uint32_t(m_x); 
    };

    // generates float from uniform distribution in [0,1) 
    __device__ __forceinline__ float randUniformFloatCO() 
    { 
        return __fdividef(__uint2float_rz(rand()), (float)0x100000000); 
    } 

    // generates float from uniform distribution in (0,1] 
    __device__ __forceinline__ float randUniformFloatOC() {return 1.0f - randUniformFloatCO();}

private:
    uint64_t m_x; // x or the current state of the rng
    uint32_t m_a; // a or the saveprime multiplier of the rng
};

/**
 * @brief allocate device memory for the mwc-rng and upload states and prime numbers
 * 
 * @param maxNumWorkitems number of prime numbers / states
 * @param MWC_RNG_x states on host
 * @param MWC_RNG_a prime numbers on host
 * @param d_MWC_RNG_x memory for states on device will be allocated here
 * @param d_MWC_RNG_a memory for primes on device will be allocated here
 */
void inline initMWCRng(int maxNumWorkitems, uint64_t* MWC_RNG_x, uint32_t* MWC_RNG_a, uint64_t** d_MWC_RNG_x,
                   uint32_t** d_MWC_RNG_a);

// function definitions
// ------------------------------------

void inline initMWCRng(int maxNumWorkitems, uint64_t* MWC_RNG_x, uint32_t* MWC_RNG_a, uint64_t** d_MWC_RNG_x,
                   uint32_t** d_MWC_RNG_a)
{
    CUDA_ERR_CHECK(cudaMalloc(d_MWC_RNG_a, maxNumWorkitems * sizeof(uint32_t)));
    CUDA_ERR_CHECK(cudaMalloc(d_MWC_RNG_x, maxNumWorkitems * sizeof(uint64_t)));

    CUDA_ERR_CHECK(cudaMemcpy(*d_MWC_RNG_a, MWC_RNG_a, maxNumWorkitems * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_ERR_CHECK(cudaMemcpy(*d_MWC_RNG_x, MWC_RNG_x, maxNumWorkitems * sizeof(uint64_t), cudaMemcpyHostToDevice));

    cudaDeviceSynchronize();
    printf("RNG is set up on CUDA gpu %d. \n", maxNumWorkitems);
}

// define a shorthand for the rng that is actually used in the code
using RngType = MultipleWithCarryRng;


#endif
