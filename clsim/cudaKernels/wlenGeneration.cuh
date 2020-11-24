/*The MIT License (MIT)

Copyright (c) 2020, Hendrik Schwanekaml, hschwanekamp@nvidia.com

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

    constexpr float wavelengthDistValues[wavelengthDistNumEntries] = {
        6.9460614465e+03f, 9.4161450667e+03f, 7.9111126641e+03f, 5.6966857194e+03f, 2.4261531760e+03f, 1.4550822538e+05f,
        5.8640971646e+05f, 1.6977055144e+06f, 3.3330250356e+06f, 5.2888747511e+06f, 6.0138643586e+06f, 6.7903489345e+06f,
        6.9089551517e+06f, 6.8157043439e+06f, 6.6069408043e+06f, 6.2897480068e+06f, 5.8715292024e+06f, 5.3891557602e+06f,
        4.9725809628e+06f, 4.6168632812e+06f, 4.2123534264e+06f, 3.7187172399e+06f, 3.2979376464e+06f, 2.9630733566e+06f,
        2.6701631226e+06f, 2.3943788910e+06f, 2.0469382067e+06f, 1.6286616843e+06f, 1.2818914213e+06f, 1.0084867047e+06f,
        8.0206182489e+05f, 6.3602259080e+05f, 5.1197464734e+05f, 4.1209066334e+05f, 3.2745735537e+05f, 2.4713577956e+05f,
        1.8606395711e+05f, 1.2297919247e+05f, 8.0123028093e+04f, 4.6665985816e+04f, 2.7703579734e+04f, 1.7255998032e+04f,
        7.4530418929e+03f,
    };

    constexpr float wavelengthDistCumulativeValues[wavelengthDistNumEntries] = {
        0.0000000000e+00f, 8.1811032566e-05f, 1.6844732122e-04f, 2.3648631314e-04f, 2.7710050761e-04f, 1.0167724004e-03f,
        4.6763621096e-03f, 1.6096938264e-02f, 4.1250591013e-02f, 8.4360089947e-02f, 1.4087378550e-01f, 2.0489485196e-01f,
        2.7339137239e-01f, 3.4201466987e-01f, 4.0912789561e-01f, 4.7361133967e-01f, 5.3441772571e-01f, 5.9072115053e-01f,
        6.4252983414e-01f, 6.9047705536e-01f, 7.3462313890e-01f, 7.7427849223e-01f, 8.0936176666e-01f, 8.4066682168e-01f,
        8.6883300407e-01f, 8.9415571414e-01f, 9.1636229963e-01f, 9.3474029908e-01f, 9.4929306461e-01f, 9.6074495524e-01f,
        9.6979769789e-01f, 9.7698811997e-01f, 9.8272810616e-01f, 9.8734843271e-01f, 9.9104617281e-01f, 9.9391913848e-01f,
        9.9608513716e-01f, 9.9763035291e-01f, 9.9864586401e-01f, 9.9927980908e-01f, 9.9965165691e-01f, 9.9987645480e-01f,
        1.0000000000e+00f,
    };

    // generates a wavelength from a random number according to the wavelength distribution
    float pregenerateWavelength(float randomNumber)
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
std::vector<float> generateWavelengthLut()
{
    std::vector<float> lut(WLEN_LUT_SIZE);
    for(int i=0; i<WLEN_LUT_SIZE; ++i)
        lut[i] = detail::pregenerateWavelength(float(i)/float(WLEN_LUT_SIZE));
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