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
    original code to generate wavelength from random numbers ported to the cpu,
    generate an entire look up table for all wlength with generateWavelengthLut()
    then use getWavelength() and pass a random number to generate a wavelength 
    lut size can be configured in the settings.cuh file
*/

#ifndef WLENGENERATION_CUH
#define WLENGENERATION_CUH

// includes
// ------------------
#include "settings.cuh"
// ------------------

namespace detail{
    // wavelength pre-generation
    constexpr int wavelengthDistNumEntries = 43;

    // generates a wavelength from a random number according to the wavelength distribution
    float pregenerateWavelength(float randomNumber, const double* wavelengthDistValues, const double* wavelengthDistCumulativeValues)
    {
        unsigned int k = 0;
        float this_acu = 0.f;  // this is 0 by definition!
        for (;;) {
            float next_acu = wavelengthDistCumulativeValues[k + 1];
            if (next_acu >= randomNumber) break;
            this_acu = next_acu;
            ++k;
        }

        // look between bins k and k+1

        const float b = wavelengthDistValues[k];
        const float x0 = float(k) * (1.0000000000e-08f) + (2.6000000000e-07f);  // _rz=="round (to) zero"
        const float slope = (wavelengthDistValues[k + 1] - b) / (1.0000000000e-08f);
        const float dy = randomNumber - this_acu;

        if ((b == 0.f) && (slope == 0.f)) {
            return x0;
        } else if (b == 0.f) {
            return x0 + sqrtf(2.f * dy / slope);

        } else if (slope == 0.f) {
            return x0 + dy / b;
        } else {
            return x0 + (sqrtf(dy * (2.f * slope) / powf(b, 2) + 1.f) - 1.f) * b / slope;
        }
    }
}

/**
 * @brief generate awavelength lookup table for use with getWavelenth()
 *          needs to be uploaded to gpu memory for use on the gpu
 *          wavelength size can be configured in settings.cuh
 */
std::vector<float> generateWavelengthLut(const double* wavelengthDistValues, const double* wavelengthDistCumulativeValues)
{
    std::vector<float> lut(WLEN_LUT_SIZE);
    for(int i=0; i<WLEN_LUT_SIZE; ++i)
        lut[i] = detail::pregenerateWavelength(float(i)/float(WLEN_LUT_SIZE), wavelengthDistValues, wavelengthDistCumulativeValues);
    return lut;
} 

/**
 * @brief get a wavelength from the look up table
 * 
 * @param randomNumber a random number in [0,1)
 * @param wlenLut pointer to the wavelength lut data (generated with generateWavelengthLut())
 * @return the generated wavelength
 */
__host__ __device__ __forceinline__ float getWavelenth(float randomNumber, const float* wlenLut)
{
    return wlenLut[int(randomNumber*WLEN_LUT_SIZE)];
}

#endif