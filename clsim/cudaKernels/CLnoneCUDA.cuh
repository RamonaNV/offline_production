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

#ifndef CLNONECUDA_CUH
#define CLNONECUDA_CUH

// mix() is a CL commonFunctions function: Returns the linear blend of x & y
// implemented as: x + (y - x) * a a must be \in  [0,1]
__device__ inline float mix(float x, float y, float a);

// clamp() is a CL commonFunctions function: Returns fmin(fmax(x, minval),
// maxval). returns x if \in [a,b], else closest boundary (a or b) Results are
// undefined if minval > maxval -> 0
__device__ inline float clamp(float x, float minval, float maxval);

// compiler says doesnt have cosnt float4 çconst float4 :
__device__ inline float4 operator*(const float4 &a, const float4 &b);

// couldnt find header for this
__device__ inline float dot(const float4 &a, const float4 &b);

__device__ inline float sign(const float4 &a);

__device__ inline float mix(float x, float y, float a) {
  return x + (y - x) * a;
}

__device__ inline float clamp(float x, float minval, float maxval) {
  return (minval <= maxval) ? fminf(fmaxf(x, minval), maxval) : 0.0;
}

__device__ inline float4 operator*(const float4 &a, const float4 &b) {
  return float4{a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w};
}

__device__ inline float dot(const float4 &a, const float4 &b) {
  return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

__device__ inline float sign(const float a) { return (a < 0.f) ? -1.0f : 1.0f; }

#endif