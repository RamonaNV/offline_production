/*The MIT License (MIT)

Copyright (c) 2020, Hendrik Schwanekamp hschwanekamp@nvidia.com, Ramona Hohl rhohl@nvidia.com

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
    collection of utility functions used in different parts of the cuda code
*/

#ifndef UTILS_CUH
#define UTILS_CUH

// includes
// ------------------
#include <cuda.h>
#include <cuda_runtime_api.h>
#include "helper_math.h" // include helper math as a utility
// ------------------

// cuda error checking

// check error in return value of cuda api call
#define CUDA_ERR_CHECK(e)              \
    if (cudaError_t(e) != cudaSuccess) \
        printf("!!! Cuda Error %s in line %d \n", cudaGetErrorString(cudaError_t(e)), __LINE__);

// chack last cuda error
#define CUDA_CHECK_LAST_ERROR()  \
    {                    \
        cudaError_t gl_err = cudaGetLastError();            \
        if (cudaError_t(gl_err) != cudaSuccess) \
            printf("!!! Cuda Error %s in line %d \n", cudaGetErrorString(cudaError_t(gl_err)), __LINE__ - 1); \
    }s
    
// math constants
constexpr float C_LIGHT = 0.299792458f;       // [m/ns]F
constexpr float RECIP_C_LIGHT = 3.33564095f;  // [ns/m]
constexpr float PI = 3.14159265359f;
constexpr float EPSILON = 0.00001f;

// math functions
__host__ __device__ __forceinline__ float sqr(const float a) { return a * a; }
__host__ __device__ __forceinline__ float mix(float x, float y, float a) { return x + (y - x) * a; }

#endif