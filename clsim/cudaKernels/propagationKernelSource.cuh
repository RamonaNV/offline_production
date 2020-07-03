
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

#ifndef PROPAGATIONKERNELSOURCE_CUH
#define PROPAGATIONKERNELSOURCE_CUH

#include <iostream>
#include <cuda.h>
#include <cuda_runtime_api.h>

//#include <dataStructCuda.cuh>
#include <opencl/mwcrng_init.h>
#include <dataStructCuda.cuh>

void finalizeCUDA();

void launch_CudaPropogate(const I3CLSimStep* __restrict__ in_steps, int nsteps,  
    uint32_t maxHitIndex, unsigned short *geoLayerToOMNumIndexPerStringSet, int ngeolayer,
    uint64_t* __restrict__  MWC_RNG_x,    uint32_t* __restrict__   MWC_RNG_a,  int sizeRNG, float& totalCudaKernelTime,const int nbenchmarks,size_t threadsPerBlock );


 
#endif  // PROPAGATIONKERNELSOURCE_CUH