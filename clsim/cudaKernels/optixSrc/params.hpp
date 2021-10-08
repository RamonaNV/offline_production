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
 
#pragma once

#include <cuda_runtime.h>
#include <optix.h>

#include <cstdint>
#include <dataStructCuda.cuh>


struct RayGenData { 
};

struct HitGroupData { 
};

struct MissData {  
};

enum GeometryType:unsigned int {DOM, CABLE, MISS};

struct Matrix
{
    float m[12];
};

struct Params {
  OptixTraversableHandle handle;
  struct I3CUDASimStep* steps;
  uint64_t* MWC_RNG_x;
  uint32_t* MWC_RNG_a;
  uint32_t nsteps, maxHitIndex;
  uint32_t* hitIndex;
  uint32_t* cableHits;
  I3CUDASimPhoton*  hitPhotons;
  float* wlenLut;
  float* zOffsetLut;
};
