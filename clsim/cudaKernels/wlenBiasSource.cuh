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

#ifndef WLENBIASSOURCE_CUH
#define WLENBIASSOURCE_CUH

//#include <CLnoneCUDA.cuh>

__device__ __forceinline__ float getWavelengthBias(float wavelength, const float* getWavelengthBias_dataShared);
__device__ __forceinline__ void getWavelengthBias_getInterpolationBinAndFraction(float wavelength, int &bin,
                                                                                 float &fraction);

__device__ float getWavelengthBias_data[43] = {
    6.3743244437e-05f, 9.3689183241e-05f, 8.5079797546e-05f, 6.6027579216e-05f, 3.0224487575e-05f, 1.9434074467e-03f,
    8.3768239340e-03f, 2.5880440148e-02f, 5.4108567428e-02f, 9.1253503761e-02f, 1.1007506732e-01f, 1.3161724587e-01f,
    1.4157858965e-01f, 1.4742781756e-01f, 1.5062904556e-01f, 1.5092873970e-01f, 1.4809603255e-01f, 1.4269860493e-01f,
    1.3806039825e-01f, 1.3425495545e-01f, 1.2815527951e-01f, 1.1824655128e-01f, 1.0949592851e-01f, 1.0262533383e-01f,
    9.6387681623e-02f, 9.0007816962e-02f, 8.0065168414e-02f, 6.6234443249e-02f, 5.4161971461e-02f, 4.4237654959e-02f,
    3.6501295149e-02f, 3.0009695504e-02f, 2.5029117564e-02f, 2.0860545296e-02f, 1.7153731691e-02f, 1.3389124157e-02f,
    1.0419306621e-02f, 7.1141064762e-03f, 4.7853569677e-03f, 2.8759740292e-03f, 1.7607982083e-03f, 1.1304887420e-03f,
    5.0300969625e-04f,
};

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