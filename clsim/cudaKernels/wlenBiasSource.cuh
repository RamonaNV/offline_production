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

#ifndef WLENBIASSOURCE_CUH
#define WLENBIASSOURCE_CUH

#include "utils.cuh"

__device__ __forceinline__ float getWavelengthBias(float wavelength, const float* getWavelengthBias_dataShared);
__device__ __forceinline__ void getWavelengthBias_getInterpolationBinAndFraction(float wavelength, int &bin,
                                                                                 float &fraction);

__device__ __forceinline__ void getWavelengthBias_getInterpolationBinAndFraction(float wavelength, int &bin,
                                                                                 float &fraction)
{
    float fbin;
    fraction = modf((wavelength - 2.6000000000e-07f) / 1.0000000000e-08f, &fbin);

    int ibin = (int)fbin;

    if ((ibin < 0) || ((ibin == 0) && (fraction < 0))) {
        ibin = 0;
        fraction = 0.f;
    } else if (ibin >= 43 - 1) {
        ibin = 43 - 2;
        fraction = 1.f;
    }

    bin = ibin;
}

__device__ __forceinline__ float getWavelengthBias(float wavelength, const float* getWavelengthBias_dataShared)
{
    int bin;
    float fraction;
    getWavelengthBias_getInterpolationBinAndFraction(wavelength, bin, fraction);

    return mix(getWavelengthBias_dataShared[bin], getWavelengthBias_dataShared[bin + 1], fraction);
}

#endif  // WLENBIASSOURCE_CUH