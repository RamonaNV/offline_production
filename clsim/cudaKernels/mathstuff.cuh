
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

#ifndef MATHSTUFF_CUH
#define MATHSTUFF_CUH

 
#include <cuda.h>
#include <cuda_runtime_api.h>

#include <math.h>
#include <math_constants.h>


typedef float real_t;
typedef float coord_t;

__host__ __device__ __forceinline__ real_t
getZenith(const coord_t &dx, const coord_t &dy, const coord_t &dz) {
  const real_t rinv = rsqrt(dx * dx + dy * dy + dz * dz);
  real_t theta = 0.;
  if (1. / rinv && abs(dz * rinv) <= 1.) {
    theta = acos(dz * rinv);
  } else {
    if (dz < 0.)
      theta = CUDART_PI_F;
  }
  if (theta < 0.)
    theta += 2. * CUDART_PI_F;
  real_t zenith = CUDART_PI_F - theta;
  if (zenith > CUDART_PI_F)
    zenith -= CUDART_PI_F - (zenith - CUDART_PI_F);

  return zenith;
}

__host__ __device__ __forceinline__ real_t
getAzimuth(const coord_t &dx, const coord_t &dy, const coord_t &dz) {
  real_t phi = 0;
  if ((dx != 0.) || (dy != 0.))
    phi = atan2(dy, dx);
  if (phi < 0.)
    phi += 2. * CUDART_PI_F;

  real_t azimuth = phi + CUDART_PI_F;
  azimuth -= (int)(azimuth / (2. * CUDART_PI_F)) * (2. * CUDART_PI_F);
  return azimuth;
}

__host__ __device__ __forceinline__ real_t
CalcTheta(const coord_t &dx, const coord_t &dy, const coord_t &dz) {

  return CUDART_PI_F - getZenith(dx, dy, dz);
}

__host__ __device__ __forceinline__ real_t CalcPhi(const coord_t &dx, const coord_t &dy,
                                          const coord_t &dz) {
  real_t phi = CUDART_PI_F + getAzimuth(dx, dy, dz);
  if (phi >= 2. * CUDART_PI_F)
    phi -= 2. * CUDART_PI_F;
  return phi;
}

#endif // MATHSTUF_CUH