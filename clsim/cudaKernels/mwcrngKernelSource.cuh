/*The MIT License (MIT)

Copyright (c) 2020, Ramona Hohl, rhohl@nvidia.com

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

#ifndef mwcrngKernelSource_CUH
#define mwcrngKernelSource_CUH

// Multiply-with-carry random number generator for OpenCL using an
// implementation along the lines of the CUDAMCML code described here:
// http://www.atomic.physics.lu.se/fileadmin/atomfysik/Biophotonics/Software/CUDAMCML.pdf

//////////////////////////////////////////////////////////////////////////////
//   Generates a random number between 0 and 1 [0,1)
//////////////////////////////////////////////////////////////////////////////
__device__ __forceinline__ float rand_MWC_co(uint64_t *x, uint32_t *a)
{
    // float rndm[10] = {0.2, 0.6, 0.43, 0.21, 0.9, 0.76, 0.1, 0.88 ,0.34, 0.78};   return rndm[*a%10];
    // Generate a random number [0,1);
    *x = (*x & 0xffffffffull) * (*a) + (*x >> 32);
    return __fdividef(__uint2float_rz((unsigned int)(*x)), (float)0x100000000);

}  // end __device__ rand_MWC_co

//////////////////////////////////////////////////////////////////////////////
//   Generates a random number between 0 and 1 (0,1]
//////////////////////////////////////////////////////////////////////////////
__device__ __forceinline__ float rand_MWC_oc(uint64_t *x, uint32_t *a) { return 1.0f - rand_MWC_co(x, a); }

#define RNG_ARGS uint64_t *rnd_x, uint32_t *rnd_a
#define RNG_ARGS_TO_CALL rnd_x, rnd_a
#define RNG_CALL_UNIFORM_CO rand_MWC_co(rnd_x, rnd_a)
#define RNG_CALL_UNIFORM_OC rand_MWC_oc(rnd_x, rnd_a)

#endif
