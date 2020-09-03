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

#include <curand_kernel.h>

struct localRngData
{
    float4 rnum;
    int numRnums{0};
};

//////////////////////////////////////////////////////////////////////////////
//   Generates a random number between 0 and 1 (0,1]
//////////////////////////////////////////////////////////////////////////////
__device__ __forceinline__ float rand_MWC_oc(curandStatePhilox4_32_10_t* thisRngState, localRngData& d)
{
    --d.numRnums;
    switch (d.numRnums)
    {
    case 2: return d.rnum.y;
    case 1: return d.rnum.z;
    case 0: return d.rnum.w;   
    case -1: 
         d.rnum = curand_uniform4(thisRngState);
         d.numRnums = 3;
         return d.rnum.x;
    }
}

//////////////////////////////////////////////////////////////////////////////
//   Generates a random number between 0 and 1 [0,1)
//////////////////////////////////////////////////////////////////////////////
__device__ __forceinline__ float rand_MWC_co(curandStatePhilox4_32_10_t* thisRngState, localRngData& d) { return 1.0f-rand_MWC_oc(thisRngState,d); }

#define RNG_ARGS curandStatePhilox4_32_10_t* thisRngState, localRngData& rngData
#define RNG_ARGS_TO_CALL thisRngState, rngData
#define RNG_CALL_UNIFORM_CO rand_MWC_co(thisRngState, rngData)
#define RNG_CALL_UNIFORM_OC rand_MWC_oc(thisRngState, rngData)

#endif
