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

#ifndef MEDIUMPROPERTYSOURCE_CUH
#define MEDIUMPROPERTYSOURCE_CUH


namespace ref{

#define STOP_PHOTONS_ON_DETECTION
#define RNG_ARGS RngType& rng
#define RNG_CALL_UNIFORM_CO rng.randUniformFloatCO()
#define RNG_CALL_UNIFORM_OC rng.randUniformFloatOC()
#define RNG_ARGS_TO_CALL rng
#define ZERO 0.0f
#define ONE 1.0f

#define MEDIUM_LAYERS 171

#define MEDIUM_MIN_WLEN 2.6500000000e-07f
#define MEDIUM_MAX_WLEN 6.7500000000e-07f

#define MEDIUM_MIN_RECIP_WLEN 1.4814814815e+06f
#define MEDIUM_MAX_RECIP_WLEN 3.7735849057e+06f

// medium layer structure:
#define MEDIUM_LAYER_BOTTOM_POS -8.5540000000e+02f
#define MEDIUM_LAYER_THICKNESS 1.0000000000e+01f

///////////////// START phase refractive index ////////////////

__device__ __forceinline__ float getPhaseRefIndex_func0(float wlen);

__device__ __forceinline__ float getPhaseRefIndex_func0(float wlen)
{
    const float n0 = 1.5574900000e+00f;
    const float n1 = -1.5798800000e+00f;
    const float n2 = 3.9999300000e+00f;
    const float n3 = -4.6827100000e+00f;
    const float n4 = 2.0935400000e+00f;

    const float x = wlen / 1e-6f;
    const float np = n0 + x * (n1 + x * (n2 + x * (n3 + x * n4)));

    return np;
}

__device__ __forceinline__ float getDispersion_func0(float wlen);

__device__ __forceinline__ float getDispersion_func0(float wlen)
{
    const float n1 = -1.5798800000e+00f;
    const float n2 = 3.9999300000e+00f;
    const float n3 = -4.6827100000e+00f;
    const float n4 = 2.0935400000e+00f;

    const float x = wlen / 1e-6f;

    const float dnp = (n1 + x * (2.f * n2 + x * (3.f * n3 + x * 4.f * n4))) / 1e-6f;

    return dnp;
}

#define FUNCTION_getPhaseRefIndex_DOES_NOT_DEPEND_ON_LAYER
__device__ __forceinline__ float getPhaseRefIndex(unsigned int layer, float wavelength);

__device__ __forceinline__ float getPhaseRefIndex(unsigned int layer, float wavelength)
{
    // phase refractive index does not have a layer structure
    return getPhaseRefIndex_func0(wavelength);
}

#define FUNCTION_getDispersion_DOES_NOT_DEPEND_ON_LAYER
__device__ __forceinline__ float getDispersion(unsigned int layer, float wavelength);

__device__ __forceinline__ float getDispersion(unsigned int layer, float wavelength)
{
    // phase refractive index does not have a layer structure
    return getDispersion_func0(wavelength);
}

///////////////// END phase refractive index ////////////////

///////////////// START group refractive index ////////////////

__device__ __forceinline__ float getGroupRefIndex_func0(float wlen);

__device__ __forceinline__ float getGroupRefIndex_func0(float wlen)
{
    const float n0 = 1.5574900000e+00f;
    const float n1 = -1.5798800000e+00f;
    const float n2 = 3.9999300000e+00f;
    const float n3 = -4.6827100000e+00f;
    const float n4 = 2.0935400000e+00f;

    const float g0 = 1.2271060000e+00f;
    const float g1 = -9.5464800000e-01f;
    const float g2 = 1.4256800000e+00f;
    const float g3 = -7.1183200000e-01f;
    const float g4 = 0.f;

    const float x = wlen / 1e-6f;
    const float np = n0 + x * (n1 + x * (n2 + x * (n3 + x * n4)));
    const float np_corr = g0 + x * (g1 + x * (g2 + x * (g3 + x * g4)));

    return np * np_corr;
}

#define FUNCTION_getGroupRefIndex_DOES_NOT_DEPEND_ON_LAYER
__device__ __forceinline__ float getGroupRefIndex(unsigned int layer, float wavelength);

__device__ __forceinline__ float getGroupRefIndex(unsigned int layer, float wavelength)
{
    // group refractive index does not have a layer structure
    return getGroupRefIndex_func0(wavelength);
}

///////////////// END group refractive index ////////////////

#ifdef FUNCTION_getGroupRefIndex_DOES_NOT_DEPEND_ON_LAYER
#define FUNCTION_getGroupVelocity_DOES_NOT_DEPEND_ON_LAYER
#endif
// group velocity from group refractive index
__device__ __forceinline__ float getGroupVelocity(unsigned int layer, float wavelength);

__device__ __forceinline__ float getGroupVelocity(unsigned int layer, float wavelength)
{
    const float c_light = 2.9979245800e-01f;
    const float n_group = getGroupRefIndex(layer, wavelength);

    return c_light / n_group;
}
///////////////// START scattering length (optimized) ////////////////

__device__ float getScatteringLength_b400[171] = {
    1.1969400000e+00f, 1.2689300000e+00f, 1.5214900000e+00f, 1.2308500000e+00f, 1.0183500000e+00f, 9.0122400000e-01f,
    5.5284400000e-01f, 2.2785500000e-01f, 1.0955600000e-01f, 8.5857400000e-02f, 1.0712800000e-01f, 1.2706000000e-01f,
    1.2067400000e-01f, 1.1404400000e-01f, 1.0687000000e-01f, 9.2694700000e-02f, 9.9412000000e-02f, 1.1814900000e-01f,
    1.4642300000e-01f, 1.8437400000e-01f, 2.2231700000e-01f, 3.1214800000e-01f, 1.9026100000e-01f, 1.8092600000e-01f,
    2.3241500000e-01f, 3.0922300000e-01f, 2.2679900000e-01f, 2.2228700000e-01f, 1.9682900000e-01f, 2.1088300000e-01f,
    1.7990500000e-01f, 1.6773000000e-01f, 1.1300900000e-01f, 1.4144100000e-01f, 9.9734000000e-02f, 1.4058400000e-01f,
    1.2786700000e-01f, 1.2522300000e-01f, 1.7784700000e-01f, 1.8057100000e-01f, 2.2776600000e-01f, 2.2600800000e-01f,
    2.9858300000e-01f, 2.8622500000e-01f, 3.8437100000e-01f, 2.9334900000e-01f, 1.8804500000e-01f, 1.7747100000e-01f,
    1.5822900000e-01f, 1.2018700000e-01f, 1.3346700000e-01f, 1.3449300000e-01f, 1.5382400000e-01f, 1.3946900000e-01f,
    1.8256700000e-01f, 1.9898400000e-01f, 1.8234200000e-01f, 2.1909000000e-01f, 2.0613900000e-01f, 3.5143100000e-01f,
    3.4280200000e-01f, 3.1711300000e-01f, 3.6292500000e-01f, 3.0966900000e-01f, 1.8837100000e-01f, 1.7623800000e-01f,
    2.3326100000e-01f, 3.3015200000e-01f, 2.8729700000e-01f, 1.8375200000e-01f, 2.0309900000e-01f, 4.3948700000e-01f,
    8.2492100000e-01f, 1.1180400000e+00f, 1.6227800000e+00f, 1.5594900000e+00f, 1.3795300000e+00f, 1.1547500000e+00f,
    8.1701200000e-01f, 9.7148800000e-01f, 1.2557500000e+00f, 1.0352200000e+00f, 4.5629700000e-01f, 5.1787200000e-01f,
    4.3512500000e-01f, 3.1947700000e-01f, 3.0481000000e-01f, 2.1828400000e-01f, 2.1999900000e-01f, 2.8487300000e-01f,
    2.4802800000e-01f, 2.5551900000e-01f, 4.7110900000e-01f, 3.7608200000e-01f, 5.2919100000e-01f, 3.8386900000e-01f,
    3.0815500000e-01f, 2.7171200000e-01f, 2.3289900000e-01f, 2.1156800000e-01f, 3.1396400000e-01f, 2.6595400000e-01f,
    2.9118700000e-01f, 3.6818700000e-01f, 4.3323700000e-01f, 4.9198900000e-01f, 4.5194600000e-01f, 4.8215300000e-01f,
    4.6678500000e-01f, 3.6282100000e-01f, 3.8418900000e-01f, 2.3666400000e-01f, 2.7425300000e-01f, 2.1962600000e-01f,
    2.7220200000e-01f, 2.5006700000e-01f, 4.0506300000e-01f, 5.9799200000e-01f, 4.7685000000e-01f, 4.3043000000e-01f,
    2.9712800000e-01f, 6.5507800000e-01f, 6.3006900000e-01f, 6.4653300000e-01f, 4.6885300000e-01f, 7.7212300000e-01f,
    6.5302300000e-01f, 5.5079600000e-01f, 2.9323400000e-01f, 2.8519500000e-01f, 2.8936800000e-01f, 2.8597700000e-01f,
    3.7204800000e-01f, 5.4050100000e-01f, 5.2094000000e-01f, 5.6446400000e-01f, 5.4535500000e-01f, 4.3253600000e-01f,
    6.0297500000e-01f, 7.2971200000e-01f, 6.0578900000e-01f, 6.5671400000e-01f, 9.6314700000e-01f, 8.1870000000e-01f,
    6.0674600000e-01f, 1.1543800000e+00f, 1.1421000000e+00f, 1.1794600000e+00f, 1.3848700000e+00f, 1.1175900000e+00f,
    9.2079500000e-01f, 9.9867200000e-01f, 1.0544100000e+00f, 1.1240100000e+00f, 1.3238500000e+00f, 1.5775900000e+00f,
    1.6638100000e+00f, 1.8036200000e+00f, 1.6263000000e+00f, 1.3691500000e+00f, 1.3407700000e+00f, 1.1310500000e+00f,
    1.1027400000e+00f, 1.2946600000e+00f, 1.7244400000e+00f, 1.7958700000e+00f, 1.6380100000e+00f, 1.5185500000e+00f,
    1.5789900000e+00f, 1.5750900000e+00f, 1.5551600000e+00f,
};

__device__ float __forceinline__ getScatteringLength(unsigned int layer, float wlen);

__device__ float __forceinline__ getScatteringLength(unsigned int layer, float wlen)
{
    const float alpha = 8.9860850573e-01f;

    return 1.f / (getScatteringLength_b400[layer] * powf(wlen * 2.5000000000e+06f, -alpha));
}

///////////////// END scattering length (optimized) ////////////////

///////////////// START absorption length (optimized) ////////////////

__device__ float getAbsorptionLength_aDust400[171] = {
    9.9900000000e+02f, 3.4163500000e-02f, 4.1300000000e-02f, 3.3092800000e-02f, 2.7146400000e-02f, 2.3892000000e-02f,
    1.4336600000e-02f, 5.6772100000e-03f, 2.6409800000e-03f, 2.0470500000e-03f, 2.5798400000e-03f, 3.0834800000e-03f,
    2.9217000000e-03f, 2.7541500000e-03f, 2.5733500000e-03f, 2.2177300000e-03f, 2.3859600000e-03f, 2.8578500000e-03f,
    3.5762100000e-03f, 4.5501800000e-03f, 5.5330800000e-03f, 7.8887000000e-03f, 5.7250700000e-03f, 3.8749600000e-03f,
    5.9264400000e-03f, 8.1376000000e-03f, 4.7020200000e-03f, 5.4081000000e-03f, 5.5004900000e-03f, 5.6618800000e-03f,
    3.9511000000e-03f, 3.3361800000e-03f, 3.2156900000e-03f, 3.9196800000e-03f, 2.6684100000e-03f, 3.5656100000e-03f,
    3.9609000000e-03f, 3.7286100000e-03f, 4.2924700000e-03f, 4.1007400000e-03f, 4.6941700000e-03f, 6.8349500000e-03f,
    5.8867200000e-03f, 9.6173000000e-03f, 1.2810300000e-02f, 1.1076200000e-02f, 4.6466300000e-03f, 4.1529300000e-03f,
    4.3853000000e-03f, 3.1557600000e-03f, 2.8060800000e-03f, 3.4751500000e-03f, 3.4406800000e-03f, 4.6345000000e-03f,
    4.6173700000e-03f, 4.9843900000e-03f, 5.1722500000e-03f, 6.7036400000e-03f, 6.3461200000e-03f, 8.6991400000e-03f,
    1.0084300000e-02f, 8.0933100000e-03f, 1.1100900000e-02f, 8.2911300000e-03f, 6.2057800000e-03f, 4.0957600000e-03f,
    6.6540400000e-03f, 6.1923400000e-03f, 6.6828300000e-03f, 4.9604200000e-03f, 6.2127400000e-03f, 1.1816700000e-02f,
    2.2866500000e-02f, 3.7382600000e-02f, 4.1686000000e-02f, 2.8713400000e-02f, 4.0505100000e-02f, 3.2661900000e-02f,
    2.0361200000e-02f, 1.9504000000e-02f, 3.9754000000e-02f, 2.2096900000e-02f, 1.6353100000e-02f, 1.0427500000e-02f,
    1.3178000000e-02f, 7.8097000000e-03f, 6.4781900000e-03f, 5.1620900000e-03f, 4.3318800000e-03f, 7.6993200000e-03f,
    5.4750700000e-03f, 6.6838700000e-03f, 1.1119600000e-02f, 1.5908900000e-02f, 9.7167400000e-03f, 1.0024200000e-02f,
    9.1608500000e-03f, 6.8498200000e-03f, 5.8824400000e-03f, 4.3388300000e-03f, 6.1388900000e-03f, 7.4258000000e-03f,
    7.7645100000e-03f, 1.0644200000e-02f, 7.7757100000e-03f, 1.0518000000e-02f, 1.9444300000e-02f, 1.5566000000e-02f,
    1.4044900000e-02f, 1.0999900000e-02f, 5.8525300000e-03f, 6.2476400000e-03f, 5.5970400000e-03f, 4.9305700000e-03f,
    6.4133800000e-03f, 8.0277700000e-03f, 8.7021700000e-03f, 1.7220700000e-02f, 1.2090000000e-02f, 1.0891000000e-02f,
    1.1674500000e-02f, 1.3377000000e-02f, 1.6312300000e-02f, 1.4357300000e-02f, 1.5135300000e-02f, 2.3892600000e-02f,
    1.7482200000e-02f, 1.4884800000e-02f, 9.9467600000e-03f, 6.8683500000e-03f, 7.2176400000e-03f, 6.2700700000e-03f,
    8.5893700000e-03f, 1.1207000000e-02f, 1.2923300000e-02f, 1.4812200000e-02f, 1.6979700000e-02f, 1.0136400000e-02f,
    2.0209100000e-02f, 2.3198700000e-02f, 2.2784400000e-02f, 1.9768900000e-02f, 1.6966800000e-02f, 2.5280000000e-02f,
    1.4658500000e-02f, 2.5593800000e-02f, 3.0603200000e-02f, 3.1650100000e-02f, 3.7432400000e-02f, 2.9917000000e-02f,
    2.4434500000e-02f, 2.6598400000e-02f, 2.8151800000e-02f, 3.0096800000e-02f, 3.5710300000e-02f, 4.2892800000e-02f,
    4.5345800000e-02f, 4.9335600000e-02f, 4.4277900000e-02f, 3.6988300000e-02f, 3.6187400000e-02f, 3.0293800000e-02f,
    2.9501800000e-02f, 3.4887900000e-02f, 4.7074200000e-02f, 4.9114100000e-02f, 4.4611300000e-02f, 4.1216800000e-02f,
    4.2932700000e-02f, 4.2821700000e-02f, 4.2255700000e-02f,
};

__device__ float getAbsorptionLength_deltaTau[171] = {
    2.7820200000e+01f,  2.7511400000e+01f,  2.7202600000e+01f,  2.6893800000e+01f,  2.6585000000e+01f,
    2.6276200000e+01f,  2.5967400000e+01f,  2.5658600000e+01f,  2.5349800000e+01f,  2.5041000000e+01f,
    2.4732200000e+01f,  2.4423400000e+01f,  2.4114600000e+01f,  2.3807100000e+01f,  2.3500700000e+01f,
    2.3195600000e+01f,  2.2891500000e+01f,  2.2588700000e+01f,  2.2286900000e+01f,  2.1986400000e+01f,
    2.1687000000e+01f,  2.1388800000e+01f,  2.1091800000e+01f,  2.0795900000e+01f,  2.0501200000e+01f,
    2.0207600000e+01f,  1.9915200000e+01f,  1.9624000000e+01f,  1.9333900000e+01f,  1.9045100000e+01f,
    1.8757300000e+01f,  1.8470800000e+01f,  1.8185400000e+01f,  1.7901100000e+01f,  1.7618000000e+01f,
    1.7336100000e+01f,  1.7055400000e+01f,  1.6775800000e+01f,  1.6497400000e+01f,  1.6220100000e+01f,
    1.5944000000e+01f,  1.5669100000e+01f,  1.5395300000e+01f,  1.5122800000e+01f,  1.4851300000e+01f,
    1.4581100000e+01f,  1.4312000000e+01f,  1.4044000000e+01f,  1.3777200000e+01f,  1.3511600000e+01f,
    1.3247200000e+01f,  1.2983900000e+01f,  1.2721800000e+01f,  1.2460800000e+01f,  1.2201100000e+01f,
    1.1942400000e+01f,  1.1685000000e+01f,  1.1428700000e+01f,  1.1173600000e+01f,  1.0919600000e+01f,
    1.0666800000e+01f,  1.0415100000e+01f,  1.0164700000e+01f,  9.9153800000e+00f,  9.6672200000e+00f,
    9.4202400000e+00f,  9.1744400000e+00f,  8.9297800000e+00f,  8.6863000000e+00f,  8.4439900000e+00f,
    8.2028200000e+00f,  7.9628300000e+00f,  7.7240000000e+00f,  7.4863400000e+00f,  7.2498500000e+00f,
    7.0145100000e+00f,  6.7803400000e+00f,  6.5473400000e+00f,  6.3154900000e+00f,  6.0848100000e+00f,
    5.8553000000e+00f,  5.6269600000e+00f,  5.3997700000e+00f,  5.1737500000e+00f,  4.9489000000e+00f,
    4.7252100000e+00f,  4.5026700000e+00f,  4.2813100000e+00f,  4.0611100000e+00f,  3.8420900000e+00f,
    3.6242200000e+00f,  3.4075200000e+00f,  3.1919700000e+00f,  2.9776000000e+00f,  2.7643900000e+00f,
    2.5523400000e+00f,  2.3414600000e+00f,  2.1317400000e+00f,  1.9231900000e+00f,  1.7158100000e+00f,
    1.5095700000e+00f,  1.3045200000e+00f,  1.1006200000e+00f,  8.9788800000e-01f,  6.9632000000e-01f,
    4.9591100000e-01f,  2.9667700000e-01f,  9.8602000000e-02f,  -9.8297000000e-02f, -2.9405200000e-01f,
    -4.8863200000e-01f, -6.8205300000e-01f, -8.7429800000e-01f, -1.0653800000e+00f, -1.2553100000e+00f,
    -1.4440600000e+00f, -1.6316500000e+00f, -1.8180800000e+00f, -2.0033600000e+00f, -2.1874500000e+00f,
    -2.3703900000e+00f, -2.5521700000e+00f, -2.7327700000e+00f, -2.9122200000e+00f, -3.0905000000e+00f,
    -3.2676100000e+00f, -3.4435600000e+00f, -3.6183500000e+00f, -3.7919600000e+00f, -3.9644200000e+00f,
    -4.1357100000e+00f, -4.3058500000e+00f, -4.4748100000e+00f, -4.6426100000e+00f, -4.8092400000e+00f,
    -4.9747200000e+00f, -5.1390200000e+00f, -5.3021600000e+00f, -5.4641400000e+00f, -5.6249500000e+00f,
    -5.7845900000e+00f, -5.9430900000e+00f, -6.1004000000e+00f, -6.2565600000e+00f, -6.4115400000e+00f,
    -6.5653700000e+00f, -6.7180300000e+00f, -6.8695400000e+00f, -7.0198700000e+00f, -7.1690400000e+00f,
    -7.3170500000e+00f, -7.4638800000e+00f, -7.6095600000e+00f, -7.7540700000e+00f, -7.8974200000e+00f,
    -8.0396000000e+00f, -8.1806200000e+00f, -8.3204700000e+00f, -8.4591500000e+00f, -8.5966800000e+00f,
    -8.7330500000e+00f, -8.8682400000e+00f, -9.0022700000e+00f, -9.1351500000e+00f, -9.2668500000e+00f,
    -9.3973800000e+00f, -9.5267500000e+00f, -9.6549700000e+00f, -9.7820100000e+00f, -9.9079000000e+00f,
    -1.0032600000e+01f,
};

__device__ __forceinline__ float getAbsorptionLength(unsigned int layer, float wlen);

__device__ __forceinline__ float getAbsorptionLength(unsigned int layer, float wlen)
{
    const float kappa = 1.0841068029e+00f;
    const float A = 6.9540903320e+03f;
    const float B = 6.6177543945e+03f;
    const float D = 6.6208071540e+02f;
    const float E = 0.f;

    const float x = wlen / 1e-9f;

    return 1.f / ((D * getAbsorptionLength_aDust400[layer] + E) * powf(x, -kappa) +
                  A * expf(-B / x) * (1.f + 0.01f * getAbsorptionLength_deltaTau[layer]));
}

///////////////// END absorption length (optimized) ////////////////

///////////////// START scattering angle distribution ////////////////

/////// BEGIN mix scattering angle generator "makeScatteringCosAngle" ////////
__device__ __forceinline__ float makeScatteringCosAngle_mix1(float rrrr__);

__device__ __forceinline__ float makeScatteringCosAngle_mix1(float rrrr__)
{
    // const float g = 9.0000000000e-01f;
    // const float beta = (1.f-g)/(1.f+g);
    const float beta = 5.2631578947e-02f;

    return clamp(2.f * powf((rrrr__), beta) - 1.f, -1.f, 1.f);
}

__device__ __forceinline__ float makeScatteringCosAngle_mix2(float rrrr__);

__device__ __forceinline__ float makeScatteringCosAngle_mix2(float rrrr__)
{
    const float g = 9.0000000000e-01f;
    const float g2 = 8.1000000000e-01f;

    // a random number [-1;+1]
    const float s = 2.f * (rrrr__)-1.f;

    const float ii = ((1.f - g2) / (1.f + g * s));
    return clamp((1.f + g2 - ii * ii) / (2.f * g), -1.f, 1.f);
}

__device__ __forceinline__ float makeScatteringCosAngle(RNG_ARGS);

__device__ __forceinline__ float makeScatteringCosAngle(RNG_ARGS)
{
    const float rr = RNG_CALL_UNIFORM_CO;
    if (rr < 3.5000000000e-01f) {
        return makeScatteringCosAngle_mix1(rr / 3.5000000000e-01f);
    } else {
        return makeScatteringCosAngle_mix2((1.f - rr) / 6.5000000000e-01f);
    }
}

///////   END mix scattering angle generator "makeScatteringCosAngle" ////////
///////////////// END scattering angle distribution ////////////////

///////////////// START directional absorption length correction function
///////////////////

__device__ __forceinline__ float getDirectionalAbsLenCorrFactor(float4 vec);

__device__ __forceinline__ float getDirectionalAbsLenCorrFactor(float4 vec)
{
    const float4 l = float4{8.5830136492e-01f, 1.0793942455e+00f, 1.0793942455e+00f, 0.f};
    const float4 rl = float4{1.1650919373e+00f, 9.2644555421e-01f, 9.2644555421e-01f, 0.f};

    const float4 n = (float4){(-6.4278760969e-01f * vec.x) + (7.6604444312e-01f * vec.y),
                              (-7.6604444312e-01f * vec.x) + (-6.4278760969e-01f * vec.y), vec.z, 0.f};
    const float4 s = n * n;

    const float nB = dot(s, rl);
    const float An = dot(s, l);

    return 2.f / ((3.0179830457e+00f - nB) * An);
}
///////////////// END directional absorption length correction function
///////////////////

///////////////// START pre-scattering direction transformation ////////////////

__device__ __forceinline__ void transformDirectionPreScatter(float4& vec);

__device__ __forceinline__ void transformDirectionPreScatter(float4& vec)
{
    vec = float4{(9.9245941798e-01f * vec.x) + (5.5392208739e-02f * vec.y) + (0.f * vec.z),
                 (5.5392208739e-02f * vec.x) + (9.7292513613e-01f * vec.y) + (0.f * vec.z),
                 (0.f * vec.x) + (0.f * vec.y) + (1.0389389999e+00f * vec.z), vec.w};

    const float norm = rsqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
    vec.x = vec.x * norm;
    vec.y = vec.y * norm;
    vec.z = vec.z * norm;
}
///////////////// END pre-scattering direction transformation ////////////////

///////////////// START post-scattering direction transformation
///////////////////

__device__ __forceinline__ void transformDirectionPostScatter(float4& vec);

__device__ __forceinline__ void transformDirectionPostScatter(float4& vec)
{
    vec = float4{(1.0108098679e+00f * vec.x) + (-5.7549125949e-02f * vec.y) + (0.f * vec.z),
                 (-5.7549125949e-02f * vec.x) + (1.0311047952e+00f * vec.y) + (0.f * vec.z),
                 (0.f * vec.x) + (0.f * vec.y) + (9.6252041756e-01f * vec.z), vec.w};

    const float norm = rsqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
    vec.x = vec.x * norm;
    vec.y = vec.y * norm;
    vec.z = vec.z * norm;
}
///////////////// END post-scattering direction transformation ////////////////

///////////////// START ice tilt z-shift ////////////////

#define getTiltZShift_data_numDistances 6
#define getTiltZShift_data_numZCoords 125
#define getTiltZShift_data_firstZCoord -5.0040000000e+02f
#define getTiltZShift_data_zCoordSpacing 1.0000000000e+01f

__device__ float getTiltZShift_data_distancesFromOriginAlongTilt[getTiltZShift_data_numDistances] = {
    -5.3141900000e+02f, -4.5488200000e+02f, 0.f, 1.6520200000e+02f, 4.4547700000e+02f, 5.2077000000e+02f,
};

__device__ float getTiltZShift_data_zCorrections[getTiltZShift_data_numDistances * getTiltZShift_data_numZCoords] = {
    4.1940800000e+01f, 4.1188100000e+01f, 4.0435400000e+01f, 3.9682700000e+01f, 3.8930000000e+01f, 3.8226800000e+01f,
    3.7612800000e+01f, 3.6835500000e+01f, 3.5834400000e+01f, 3.4466000000e+01f, 3.3040500000e+01f, 3.1387800000e+01f,
    3.0063700000e+01f, 2.8294400000e+01f, 2.6889800000e+01f, 2.5276500000e+01f, 2.3713000000e+01f, 2.2398500000e+01f,
    2.1300300000e+01f, 2.0489600000e+01f, 1.9996700000e+01f, 1.9668800000e+01f, 1.9248100000e+01f, 1.8442100000e+01f,
    1.7100100000e+01f, 1.5172000000e+01f, 1.2801600000e+01f, 1.0474600000e+01f, 8.5759800000e+00f, 7.2873000000e+00f,
    6.2278000000e+00f, 5.3461300000e+00f, 4.5783900000e+00f, 3.8055300000e+00f, 3.0093500000e+00f, 2.3087200000e+00f,
    1.4637100000e+00f, 2.0222200000e-01f, -1.0952900000e+00f, -2.5619500000e+00f, -3.9427300000e+00f,
    -5.0178700000e+00f, -5.4536700000e+00f, -5.3571300000e+00f, -5.2660400000e+00f, -5.3108200000e+00f,
    -5.2127200000e+00f, -4.9138100000e+00f, -4.4959600000e+00f, -4.1450000000e+00f, -3.8552900000e+00f,
    -3.7020000000e+00f, -3.7955100000e+00f, -4.0782500000e+00f, -4.6558700000e+00f, -5.2153600000e+00f,
    -5.5854600000e+00f, -5.5321400000e+00f, -4.9599100000e+00f, -3.8815000000e+00f, -2.7985700000e+00f,
    -1.8746300000e+00f, -1.3046200000e+00f, -1.1033300000e+00f, -1.3215500000e+00f, -1.8377400000e+00f,
    -2.6842900000e+00f, -3.7756200000e+00f, -4.7011800000e+00f, -5.0942400000e+00f, -4.9680400000e+00f,
    -4.9205100000e+00f, -5.1994700000e+00f, -5.8033300000e+00f, -6.2596900000e+00f, -6.3922200000e+00f,
    -6.3200000000e+00f, -6.1547600000e+00f, -5.8214600000e+00f, -5.4075000000e+00f, -5.1543100000e+00f,
    -5.1180000000e+00f, -5.0170000000e+00f, -4.9190000000e+00f, -5.0485700000e+00f, -5.2481300000e+00f,
    -5.7312900000e+00f, -6.7902200000e+00f, -8.2558800000e+00f, -1.0060500000e+01f, -1.1956900000e+01f,
    -1.3333400000e+01f, -1.3765000000e+01f, -1.3559400000e+01f, -1.3049800000e+01f, -1.2216300000e+01f,
    -1.1407600000e+01f, -1.0610700000e+01f, -9.9209400000e+00f, -9.3718900000e+00f, -8.8414300000e+00f,
    -8.3566700000e+00f, -7.8842900000e+00f, -7.4748500000e+00f, -7.1267300000e+00f, -6.7642300000e+00f,
    -6.3565400000e+00f, -5.8549100000e+00f, -5.3322600000e+00f, -4.7218900000e+00f, -4.1681000000e+00f,
    -3.7109500000e+00f, -3.3507700000e+00f, -3.0282500000e+00f, -2.7884500000e+00f, -2.5435300000e+00f,
    -2.3059200000e+00f, -2.0915700000e+00f, -1.7749000000e+00f, -1.4882700000e+00f, -1.1612600000e+00f,
    -9.1411800000e-01f, -7.0663400000e-01f, -6.8212100000e-01f, -7.7505100000e-01f,
    // distances[0]
    3.9415100000e+01f, 3.8776800000e+01f, 3.8138500000e+01f, 3.7500200000e+01f, 3.6861900000e+01f, 3.6356800000e+01f,
    3.5825700000e+01f, 3.5064100000e+01f, 3.4135400000e+01f, 3.3062200000e+01f, 3.1928900000e+01f, 3.0467900000e+01f,
    2.9030000000e+01f, 2.7913500000e+01f, 2.6728900000e+01f, 2.5398200000e+01f, 2.4027800000e+01f, 2.2756700000e+01f,
    2.1611300000e+01f, 2.0789600000e+01f, 2.0236200000e+01f, 1.9886700000e+01f, 1.9620600000e+01f, 1.9141700000e+01f,
    1.8326600000e+01f, 1.6818200000e+01f, 1.4900200000e+01f, 1.2918100000e+01f, 1.1178300000e+01f, 9.8811100000e+00f,
    8.8585700000e+00f, 8.0278300000e+00f, 7.2214000000e+00f, 6.4912900000e+00f, 5.8097900000e+00f, 5.0789400000e+00f,
    4.3429000000e+00f, 3.4948400000e+00f, 2.4755600000e+00f, 1.2547200000e+00f, 4.1111100000e-02f, -1.0765900000e+00f,
    -1.7046900000e+00f, -1.7111800000e+00f, -1.5258800000e+00f, -1.3452500000e+00f, -1.2817600000e+00f,
    -1.1294100000e+00f, -9.3372500000e-01f, -7.4843100000e-01f, -5.3274500000e-01f, -3.9323200000e-01f,
    -4.5877600000e-01f, -7.2416700000e-01f, -1.1963200000e+00f, -1.7763200000e+00f, -2.2606200000e+00f,
    -2.4150000000e+00f, -2.2176200000e+00f, -1.4636400000e+00f, -5.7363600000e-01f, 2.5222200000e-01f,
    8.8428600000e-01f, 1.2012900000e+00f, 1.1382500000e+00f, 7.5000000000e-01f, 9.8817200000e-02f, -8.2054900000e-01f,
    -1.7420400000e+00f, -2.2720600000e+00f, -2.3550000000e+00f, -2.2848500000e+00f, -2.4189800000e+00f,
    -2.8252100000e+00f, -3.1174200000e+00f, -3.2851500000e+00f, -3.3250000000e+00f, -3.1886300000e+00f,
    -3.0068900000e+00f, -2.6572500000e+00f, -2.5528300000e+00f, -2.6709700000e+00f, -2.6086100000e+00f,
    -2.5210000000e+00f, -2.5490000000e+00f, -2.6053500000e+00f, -2.8826300000e+00f, -3.6145700000e+00f,
    -4.7631000000e+00f, -6.3781400000e+00f, -8.0304700000e+00f, -9.3966700000e+00f, -1.0182200000e+01f,
    -1.0184600000e+01f, -9.7360400000e+00f, -9.1120600000e+00f, -8.4541100000e+00f, -7.8232700000e+00f,
    -7.2681000000e+00f, -6.7795200000e+00f, -6.3338100000e+00f, -5.9020800000e+00f, -5.5068900000e+00f,
    -5.1670900000e+00f, -4.8777700000e+00f, -4.5700000000e+00f, -4.2758300000e+00f, -3.9238500000e+00f,
    -3.5176200000e+00f, -3.0930800000e+00f, -2.7055800000e+00f, -2.3797100000e+00f, -2.0847100000e+00f,
    -1.8010700000e+00f, -1.5856900000e+00f, -1.3389300000e+00f, -1.1520000000e+00f, -1.0020400000e+00f,
    -7.9549000000e-01f, -5.7294100000e-01f, -3.4058800000e-01f, -2.1411800000e-01f, -7.4000000000e-02f,
    -5.6868700000e-02f, -2.3632700000e-01f,
    // distances[1]
    0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
    0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
    0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
    0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
    0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
    0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
    // distances[2]
    -1.8408100000e+01f, -1.7918600000e+01f, -1.7458300000e+01f, -1.7018100000e+01f, -1.6292400000e+01f,
    -1.5675600000e+01f, -1.4943400000e+01f, -1.4316100000e+01f, -1.3746100000e+01f, -1.3064300000e+01f,
    -1.2610000000e+01f, -1.2022800000e+01f, -1.1366300000e+01f, -1.0709800000e+01f, -9.9580700000e+00f,
    -9.1067000000e+00f, -8.2883500000e+00f, -7.4922200000e+00f, -6.7316800000e+00f, -6.1067900000e+00f,
    -5.6176200000e+00f, -5.2151900000e+00f, -4.9111800000e+00f, -4.6321400000e+00f, -4.2978800000e+00f,
    -3.8890500000e+00f, -3.4414300000e+00f, -2.9538100000e+00f, -2.5065400000e+00f, -2.1263100000e+00f,
    -1.7940400000e+00f, -1.4382700000e+00f, -1.0804800000e+00f, -6.9596200000e-01f, -2.9884600000e-01f,
    7.7619000000e-02f, 3.9601900000e-01f, 7.1269200000e-01f, 1.0357700000e+00f, 1.3484500000e+00f, 1.6973100000e+00f,
    2.1300000000e+00f, 2.5338500000e+00f, 2.7670000000e+00f, 2.8080000000e+00f, 2.6718400000e+00f, 2.4238800000e+00f,
    2.1667300000e+00f, 1.9936400000e+00f, 1.9148500000e+00f, 1.8754500000e+00f, 1.8680000000e+00f, 1.9191100000e+00f,
    2.0606900000e+00f, 2.1547500000e+00f, 2.3790200000e+00f, 2.6054900000e+00f, 2.7844600000e+00f, 2.8705900000e+00f,
    2.8835400000e+00f, 2.6784500000e+00f, 2.3112500000e+00f, 1.9070800000e+00f, 1.5805200000e+00f, 1.3724200000e+00f,
    1.2890000000e+00f, 1.3468300000e+00f, 1.5193200000e+00f, 1.7668900000e+00f, 2.0999000000e+00f, 2.2996100000e+00f,
    2.3560000000e+00f, 2.3340000000e+00f, 2.4181200000e+00f, 2.5488100000e+00f, 2.6003000000e+00f, 2.5800000000e+00f,
    2.5047500000e+00f, 2.4744400000e+00f, 2.2851000000e+00f, 2.0585700000e+00f, 1.9676200000e+00f, 2.0535300000e+00f,
    2.2005900000e+00f, 2.3686100000e+00f, 2.4950000000e+00f, 2.5140000000e+00f, 2.5108100000e+00f, 2.6270600000e+00f,
    2.9938100000e+00f, 3.5720600000e+00f, 4.1790600000e+00f, 4.7175000000e+00f, 5.1280600000e+00f, 5.4035300000e+00f,
    5.5201000000e+00f, 5.5050000000e+00f, 5.4200000000e+00f, 5.4552500000e+00f, 5.3855600000e+00f, 5.3209100000e+00f,
    5.3017200000e+00f, 5.1789800000e+00f, 5.0055100000e+00f, 4.8565300000e+00f, 4.6463300000e+00f, 4.4116300000e+00f,
    4.1846400000e+00f, 3.9062900000e+00f, 3.6126500000e+00f, 3.3196900000e+00f, 3.1055100000e+00f, 2.8926300000e+00f,
    2.7534700000e+00f, 2.6067700000e+00f, 2.4228600000e+00f, 2.2542400000e+00f, 2.0830600000e+00f, 1.9441400000e+00f,
    1.7524500000e+00f, 1.6077800000e+00f, 1.5440000000e+00f, 1.5480000000e+00f, 1.5027300000e+00f, 1.6427500000e+00f,
    // distances[3]
    -1.3245500000e+01f, -1.3274000000e+01f, -1.3511700000e+01f, -1.3923600000e+01f, -1.4157900000e+01f,
    -1.3961400000e+01f, -1.3306700000e+01f, -1.2520900000e+01f, -1.1525000000e+01f, -1.0513200000e+01f,
    -9.5609100000e+00f, -8.7030300000e+00f, -7.8525700000e+00f, -6.9609100000e+00f, -6.0798200000e+00f,
    -5.0436400000e+00f, -3.9395700000e+00f, -2.7030400000e+00f, -1.4542100000e+00f, -2.5660700000e-01f,
    6.8818200000e-01f, 1.4083000000e+00f, 1.8775700000e+00f, 2.1652900000e+00f, 2.2650000000e+00f, 2.3587100000e+00f,
    2.5103900000e+00f, 2.8614300000e+00f, 3.3195200000e+00f, 3.7766700000e+00f, 4.1742300000e+00f, 4.5559600000e+00f,
    5.0421500000e+00f, 5.6545300000e+00f, 6.1557100000e+00f, 6.7073600000e+00f, 7.2338100000e+00f, 7.6530800000e+00f,
    7.9092100000e+00f, 8.1569200000e+00f, 8.5396200000e+00f, 9.0113100000e+00f, 9.7300000000e+00f, 1.0405200000e+01f,
    1.0717900000e+01f, 1.0656000000e+01f, 1.0425800000e+01f, 9.8736200000e+00f, 9.2619100000e+00f, 8.7215800000e+00f,
    8.2588700000e+00f, 7.9626500000e+00f, 7.7653500000e+00f, 7.7399000000e+00f, 7.7339600000e+00f, 7.9876900000e+00f,
    8.3713500000e+00f, 8.7557100000e+00f, 9.0940800000e+00f, 9.2860000000e+00f, 9.2310300000e+00f, 8.6525800000e+00f,
    7.7607700000e+00f, 6.7637000000e+00f, 5.8816100000e+00f, 5.3196900000e+00f, 5.1570000000e+00f, 5.3173800000e+00f,
    5.6742300000e+00f, 6.2010300000e+00f, 6.9064200000e+00f, 7.2996100000e+00f, 7.3460000000e+00f, 7.2310100000e+00f,
    7.1820400000e+00f, 7.0932700000e+00f, 6.9228600000e+00f, 6.6588700000e+00f, 6.4372200000e+00f, 6.1454600000e+00f,
    5.8708200000e+00f, 5.7470000000e+00f, 5.9474800000e+00f, 6.0874300000e+00f, 6.3416500000e+00f, 6.5060000000e+00f,
    6.4698000000e+00f, 6.3087900000e+00f, 6.3399000000e+00f, 6.8896300000e+00f, 7.9291200000e+00f, 9.2317400000e+00f,
    1.0434400000e+01f, 1.1522900000e+01f, 1.2317900000e+01f, 1.2727100000e+01f, 1.2902000000e+01f, 1.2814700000e+01f,
    1.2694600000e+01f, 1.2535100000e+01f, 1.2389600000e+01f, 1.2181000000e+01f, 1.1899800000e+01f, 1.1525800000e+01f,
    1.1105800000e+01f, 1.0658400000e+01f, 1.0219600000e+01f, 9.7984200000e+00f, 9.3883300000e+00f, 9.0330900000e+00f,
    8.6622900000e+00f, 8.3403100000e+00f, 8.0473500000e+00f, 7.7695800000e+00f, 7.3810400000e+00f, 7.0093800000e+00f,
    6.6763900000e+00f, 6.2925000000e+00f, 6.0104100000e+00f, 5.6363200000e+00f, 5.3093800000e+00f, 5.1116300000e+00f,
    5.0000000000e+00f, 5.0349500000e+00f, 5.0258800000e+00f,
    // distances[4]
    -1.2580000000e+00f, -1.2224800000e+00f, -1.3575000000e+00f, -2.0852200000e+00f, -2.9345200000e+00f,
    -3.5829000000e+00f, -3.9450000000e+00f, -3.9610900000e+00f, -3.3956900000e+00f, -2.3941100000e+00f,
    -1.3351800000e+00f, -4.2740700000e-01f, 3.7587200000e-01f, 1.3146800000e+00f, 2.2790900000e+00f, 3.2083800000e+00f,
    4.1636400000e+00f, 5.2246900000e+00f, 6.4326100000e+00f, 7.7291300000e+00f, 8.8972600000e+00f, 9.9436400000e+00f,
    1.0626200000e+01f, 1.1062400000e+01f, 1.1096700000e+01f, 1.0891200000e+01f, 1.0665700000e+01f, 1.0558300000e+01f,
    1.0550800000e+01f, 1.0708200000e+01f, 1.0925100000e+01f, 1.1225100000e+01f, 1.1587500000e+01f, 1.1968500000e+01f,
    1.2453800000e+01f, 1.2965800000e+01f, 1.3486700000e+01f, 1.3958600000e+01f, 1.4331900000e+01f, 1.4524100000e+01f,
    1.4620200000e+01f, 1.5039600000e+01f, 1.5550400000e+01f, 1.6360600000e+01f, 1.7005200000e+01f, 1.7248800000e+01f,
    1.7182000000e+01f, 1.6806300000e+01f, 1.6126800000e+01f, 1.5259300000e+01f, 1.4467600000e+01f, 1.3778900000e+01f,
    1.3349800000e+01f, 1.3094600000e+01f, 1.3011200000e+01f, 1.3116400000e+01f, 1.3418600000e+01f, 1.3713300000e+01f,
    1.3964600000e+01f, 1.4272200000e+01f, 1.4318800000e+01f, 1.3920300000e+01f, 1.2961500000e+01f, 1.1684500000e+01f,
    1.0436600000e+01f, 9.5161700000e+00f, 9.0430000000e+00f, 9.0290000000e+00f, 9.2843700000e+00f, 9.6907500000e+00f,
    1.0309400000e+01f, 1.0864600000e+01f, 1.1020900000e+01f, 1.0896700000e+01f, 1.0808000000e+01f, 1.0764000000e+01f,
    1.0516600000e+01f, 1.0353700000e+01f, 1.0035200000e+01f, 9.5310400000e+00f, 9.0005300000e+00f, 8.5453100000e+00f,
    8.5230000000e+00f, 8.8215100000e+00f, 9.1569200000e+00f, 9.5202900000e+00f, 9.6158600000e+00f, 9.3990700000e+00f,
    9.2560000000e+00f, 9.6337000000e+00f, 1.0583500000e+01f, 1.1936800000e+01f, 1.3331800000e+01f, 1.4635300000e+01f,
    1.5674000000e+01f, 1.6330000000e+01f, 1.6652500000e+01f, 1.6737000000e+01f, 1.6689600000e+01f, 1.6631000000e+01f,
    1.6576500000e+01f, 1.6443300000e+01f, 1.6202200000e+01f, 1.5777900000e+01f, 1.5318500000e+01f, 1.4761600000e+01f,
    1.4168900000e+01f, 1.3564000000e+01f, 1.2956300000e+01f, 1.2360500000e+01f, 1.1782600000e+01f, 1.1260500000e+01f,
    1.0814400000e+01f, 1.0421700000e+01f, 9.9258800000e+00f, 9.4563200000e+00f, 8.9841700000e+00f, 8.4887600000e+00f,
    8.0100000000e+00f, 7.5657900000e+00f, 7.1154200000e+00f, 6.8491900000e+00f, 6.7062900000e+00f, 6.6785100000e+00f,
    6.6570000000e+00f,
    // distances[5]
};

__device__ __forceinline__ float getTiltZShift(const float4& vec);

__device__ __forceinline__ float getTiltZShift(const float4& vec)
{
    const float z_rescaled = (vec.z - getTiltZShift_data_firstZCoord) / getTiltZShift_data_zCoordSpacing;
    const int k = min(max(__float2int_rd(z_rescaled), 0), getTiltZShift_data_numZCoords - 2);

    const float fraction_z_above = z_rescaled - float(k);
    const float fraction_z_below = 1. - fraction_z_above;

    const float nr = -7.0710678119e-01f * vec.x + -7.0710678119e-01f * vec.y;

    for (int j = 1; j < getTiltZShift_data_numDistances; j++) {
        const float thisDist = getTiltZShift_data_distancesFromOriginAlongTilt[j];
        if ((nr < thisDist) || (j == getTiltZShift_data_numDistances - 1)) {
            const float previousDist = getTiltZShift_data_distancesFromOriginAlongTilt[j - 1];
            const float thisDistanceBinWidth = thisDist - previousDist;

            const float frac_at_lower = (thisDist - nr) / thisDistanceBinWidth;
            const float frac_at_upper = 1. - frac_at_lower;

            const float val_at_lower =
                (getTiltZShift_data_zCorrections[(j - 1) * getTiltZShift_data_numZCoords + k + 1] * fraction_z_above +
                 getTiltZShift_data_zCorrections[(j - 1) * getTiltZShift_data_numZCoords + k] * fraction_z_below);
            const float val_at_upper =
                (getTiltZShift_data_zCorrections[j * getTiltZShift_data_numZCoords + k + 1] * fraction_z_above +
                 getTiltZShift_data_zCorrections[j * getTiltZShift_data_numZCoords + k] * fraction_z_below);

            return (val_at_upper * frac_at_upper + val_at_lower * frac_at_lower);
        }
    }
    return 0;
}
///////////////// END ice tilt z-shift ////////////////

}
#endif  // MEDIUMPROPERTYSOURCE_CUH