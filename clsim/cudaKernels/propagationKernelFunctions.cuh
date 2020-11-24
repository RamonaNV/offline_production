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
    contains all helper functions used to actually run the simulation
*/

#ifndef PROPAGATIONKERNELFUNCTIUONS_CUH
#define PROPAGATIONKERNELFUNCTIUONS_CUH

// includes
// ------------------
#include <cuda.h>
#include <cuda_runtime_api.h>

#include "settings.cuh"
#include "dataStructCuda.cuh"
#include "utils.cuh"
#include "rng.cuh"
#include "domGeoData.cuh"
#include "wlenBiasSource.cuh"
#include "zOffsetHandling.cuh"
#include "wlenGeneration.cuh"
// ------------------

// calculates the step direction
__device__ __forceinline__ float3 calculateStepDir(const I3CLSimStepCuda& step);

/**
 * @brief Creates a single photon to be propagated
 * @param step the step to create the photon from
 * @param stepDir step direction to create the photon ( calculated using calculateStepDir() )
 * @param wlenLut look up table to generate wavelength from random numbers
 * @param rng the rng to use for generating this photon, 4 rng values are computed
 */
__device__ __forceinline__ I3CLInitialPhoton createPhoton(const I3CLSimStepCuda& step, const float3& stepDir, const float* wlenLut, RngType& rng);

/**
 * @brief  propgates a single photon once
 * @param ph the photon to propagate
 * @param distancePropagated the distance the photon was propagated during this iteration
 * @param rng the rng to use for generating this photon, 1 rng values is computed
 * @param scatteringLength scattering length look up table, which can be in global or shared memory
 * @param absorptionDust dust absorption look up table, which can be in global or shared memory
 * @param absorptionTauDelta absorption Tau-Delta look up table, which can be in global or shared memory
 * @param zOffsetLut lut containing zOffset values 
 * @return true if the photon was absorbed
 */
__device__ __forceinline__ bool propPhoton(I3CLPhoton& ph, float& distancePropagated, RngType rng, const float* scatteringLength, const loat* absorptionDust, const float* absorptionTauDelta, const float* zOffsetLut);

/**
 * @brief moves a photon along its track by the propagated distance
 * @param ph the photon to move
 * @param distancePropagated the distance the photon was propagated this iteration
 */
__device__ __forceinline__  void updatePhotonTrack(I3CLPhoton& ph, float distancePropagated);

/**
 * @brief scatters a photon
 * @param ph the photon to scatter
 * @param rng the random number generator 
 */
__device__ __forceinline__  void scatterPhoton(I3CLPhoton& ph, RngType& rng);

// helper functions and subfunctions
// --------------------------------------------------------------------------------------
namespace detail {
    // helper functions to compute wavelength dependent properties
    // -------------------

    // for second argument "x" pass wlen*1e6f 
    __device__ __forceinline__ float getPhaseRefIndex(float wlen, float x)
    {
        constexpr float n0 = 1.5574900000e+00f;
        constexpr float n1 = -1.5798800000e+00f;
        constexpr float n2 = 3.9999300000e+00f;
        constexpr float n3 = -4.6827100000e+00f;
        constexpr float n4 = 2.0935400000e+00f;

        const float np = n0 + x * (n1 + x * (n2 + x * (n3 + x * n4)));
        return np;
    }

    // use this if "x" is not known already
    __device__ __forceinline__ float getPhaseRefIndex(float wlen)
    {
        return getPhaseRefIndex(wlen, wlen*1e6f);
    }

    __device__ __forceinline__ float getDispersion(float wlen)
    {
        constexpr float n1 = -1.5798800000e+00f;
        constexpr float n2 = 3.9999300000e+00f;
        constexpr float n3 = -4.6827100000e+00f;
        constexpr float n4 = 2.0935400000e+00f;

        const float x = wlen * 1e6f;
        const float dnp = (n1 + x * (2.f * n2 + x * (3.f * n3 + x * 4.f * n4))) * 1e6f;

        return dnp;
    }

    // for the second argument, pass the result of getPhaseRefIndex() to avoid recalculation,
    // for "x" pass wlen*1e6f
    __device__ __forceinline__ float getGroupRefIndex(float wlen, float phaseRefIndex, float x)
    {
        constexpr float g0 = 1.2271060000e+00f;
        constexpr float g1 = -9.5464800000e-01f;
        constexpr float g2 = 1.4256800000e+00f;
        constexpr float g3 = -7.1183200000e-01f;
        constexpr float g4 = 0.f;

        const float np_corr = g0 + x * (g1 + x * (g2 + x * (g3 + x * g4)));

        return phaseRefIndex * np_corr;
    }

    // use if phaseRefIndex and x are not known already
    __device__ __forceinline__ float getGroupRefIndex(float wlen)
    {
        const float x = wlen*1e6f;
        return getGroupRefIndex(wlen,getPhaseRefIndex(wlen,x),x);
    }

    // scattering helper functions
    // -------------------

    // scattering by rotating the direction vector by a random amount
    __device__ __forceinline__ void scatterDirectionByAngle(float cosa, float sina, float3 &direction, float randomNumber)
    {
        // randomize direction of scattering (rotation around old direction axis)
        const float b = 2.0f * PI * randomNumber;

        // Rotate new direction into absolute frame of reference
        const float cosb = cosf(b);
        const float sinb = sinf(b);
        const float sinth = sqrtf(max(0.0f, 1.0f - sqr(direction.z)));

        if (sinth > 0.f) {  // Current direction not vertical, so rotate
            const float3 oldDir = direction;
            direction.x = oldDir.x * cosa - ((oldDir.y * cosb + oldDir.z * oldDir.x * sinb) * sina) / sinth;
            direction.y = oldDir.y * cosa + ((oldDir.x * cosb - oldDir.z * oldDir.y * sinb) * sina) / sinth;
            direction.z = oldDir.z * cosa + sina * sinb * sinth;
        } else {  // Current direction is vertical, so this is trivial
            direction.x = sina * cosb;
            direction.y = sina * sinb;
            direction.z = cosa * copysignf(1.0f, direction.z);
        }

        direction = normalize(direction);
    }

    // perform pre-scatter transformation on photon direction
    __device__ __forceinline__ void transformDirectionPreScatter(float3& dir)
    {
        dir = float3{(9.9245941798e-01f * dir.x) + (5.5392208739e-02f * dir.y) + (0.f * dir.z),
                    (5.5392208739e-02f * dir.x) + (9.7292513613e-01f * dir.y) + (0.f * dir.z),
                    (0.f * dir.x) + (0.f * dir.y) + (1.0389389999e+00f * dir.z)};
        dir = normalize(dir);
    }

    // perform post-scatter transformation on photon direction
    __device__ __forceinline__ void transformDirectionPostScatter(float3& dir)
    {
        dir = float4{(1.0108098679e+00f * dir.x) + (-5.7549125949e-02f * dir.y) + (0.f * dir.z),
                    (-5.7549125949e-02f * dir.x) + (1.0311047952e+00f * dir.y) + (0.f * dir.z),
                    (0.f * dir.x) + (0.f * dir.y) + (9.6252041756e-01f * dir.z)};

            dir = normalize(dir);
    }

    // compute first variant of the scattering cosine
    __device__ __forceinline__ float makeScatteringCosAngle_mix1(float rrrr__)
    {
        // const float g = 9.0000000000e-01f;
        // const float beta = (1.f-g)/(1.f+g);
        const float beta = 5.2631578947e-02f;

        return clamp(2.f * powf((rrrr__), beta) - 1.f, -1.f, 1.f);
    }

    // compute second variant of the scattering cosine
    __device__ __forceinline__ float makeScatteringCosAngle_mix2(float rrrr__)
    {
        const float g = 9.0000000000e-01f;
        const float g2 = 8.1000000000e-01f;

        // a random number [-1;+1]
        const float s = 2.f * (rrrr__)-1.f;

        const float ii = ((1.f - g2) / (1.f + g * s));
        return clamp((1.f + g2 - ii * ii) / (2.f * g), -1.f, 1.f);
    }

    // compute the scattering cosine by selecting between variant 1 and 2
    __device__ __forceinline__ float makeScatteringCosAngle(float randomNumberCO)
    {
        if (randomNumberCO < 3.5000000000e-01f) {
            return makeScatteringCosAngle_mix1(randomNumberCO / 3.5000000000e-01f);
        } else {
            return makeScatteringCosAngle_mix2((1.f - randomNumberCO) / 6.5000000000e-01f);
        }
    }

    // propagation helper funtions
    // -------------------

    // compute the absorbtion length correction factor
    __device__ __forceinline__ float getDirectionalAbsLenCorrFactor(float3 vec)
    {
        const float3 l = float3{8.5830136492e-01f, 1.0793942455e+00f, 1.0793942455e+00f};
        const float3 rl = float3{1.1650919373e+00f, 9.2644555421e-01f, 9.2644555421e-01f};

        const float3 n = float3{(-6.4278760969e-01f * vec.x) + (7.6604444312e-01f * vec.y),
                                (-7.6604444312e-01f * vec.x) + (-6.4278760969e-01f * vec.y), vec.z};
        const float3 s = n * n;

        const float nB = dot(s, rl);
        const float An = dot(s, l);

        return 2.f / ((3.0179830457e+00f - nB) * An);
    }

    // compute scattering length in current layer from lut values
    __device__ float __forceinline__ getScatteringLength(unsigned int layer, float wlen, const float* scatteringLengthLut)
    {
        const float alpha = 8.9860850573e-01f;
        return 1.f / (scatteringLengthLut[layer] * powf(wlen * 2.5000000000e+06f, -alpha));
    }

    // compute absorption length in current layer from lut values
    __device__ __forceinline__ float getAbsorptionLength(unsigned int layer, float wlen, const float* absorptionADustLut, const float* absorptionDeltaTauLut)
    {
        const float kappa = 1.0841068029e+00f;
        const float A = 6.9540903320e+03f;
        const float B = 6.6177543945e+03f;
        const float D = 6.6208071540e+02f;
        const float E = 0.f;

        const float x = wlen / 1e-9f;

        return 1.f / ((D * absorptionADustLut[layer] + E) * powf(x, -kappa) +
                    A * expf(-B / x) * (1.f + 0.01f * absorptionDeltaTauLut[layer]));
    }
}

// function definitions for the major functions
// --------------------------------------------------------------------------------------

__device__ __forceinline__ float3 calculateStepDir(const I3CLSimStepCuda& step)
{
        const float rho = sinf(step.dirAndLengthAndBeta.x);       // sin(theta)
        return float3{rho * cosf(step.dirAndLengthAndBeta.y),  // rho*cos(phi)
                         rho * sinf(step.dirAndLengthAndBeta.y),  // rho*sin(phi)
                         cosf(step.dirAndLengthAndBeta.x)};        // cos(phi)
}

__device__ __forceinline__ I3CLInitialPhoton createPhoton(const I3CLSimStepCuda& step, const float3& stepDir, const float* wlenLut, RngType& rng)
{
    // float4 randomNumbers = float4{0.12454854,0.99568,0.4877858,0.24784564};
    float4 randomNumbers = {rng.randUniformFloatCO(), rng.randUniformFloatCO(), rng.randUniformFloatCO(), rng.randUniformFloatOC()};

    I3CLInitialPhoton ph;

    // move along the step direction a random amount
    const float shiftMultiplied = step.dirAndLengthAndBeta.z * randomNumbers.x;
    const float inverseParticleSpeed = 1.f / (C_LIGHT * step.dirAndLengthAndBeta.w);
    ph.pos = float3{step.posAndTime.x + stepDir.x * shiftMultiplied, step.posAndTime.y + stepDir.y * shiftMultiplied,
                    step.posAndTime.z + stepDir.z * shiftMultiplied};
    ph.time = step.posAndTime.w + inverseParticleSpeed * shiftMultiplied;

    // generate a wavelength
    ph.wlen = getWlen(randomNumbers.y,wlenLut);

    // calculate phase and group ref index
    const float x = ph.wlen * 1e6;
    const float phaseRefIndex = getPhaseRefIndex(ph.wlen,x);
    const float groupRefIndex = getGroupRefIndex(ph.wlen,phaseRefIndex,x);

    const float cosCherenkov = min( 1.0f, 1.0f / (step.dirAndLengthAndBeta.w * phaseRefIndex));  // cos theta = 1/(beta*n)
    const float sinCherenkov = sqrtf(1.0f - sqr(cosCherenkov));
    
    // determine the photon direction
    // start with the track direction and rotate to cherenkov emission direction
    ph.dir = stepDir;
    scatterDirectionByAngle(cosCherenkov, sinCherenkov, ph.dir, randomNumbers.z);
    
    // calc inverse group velocity and set an initial absorption length
    ph.invGroupvel = groupRefIndex * RECIP_C_LIGHT; // refIndex * (1/c_light) <=> 1 / (c_light / refIndex)
    ph.absLength = -logf(randomNumbers.w);

    return ph;
}

__device__ __forceinline__ bool propPhoton(I3CLPhoton& ph, float& distancePropagated, RngType rng, const float* scatteringLength, const loat* absorptionDust, const float* absorptionTauDelta, const float* zOffsetLut)
{ 
    const float effective_z = ph.pos.z - getZOffset(ph.pos, zOffsetLut);
    const int currentPhotonLayer = min(max(findLayerForGivenZPos(effective_z), 0), MEDIUM_LAYERS - 1);
    const float photon_dz = ph.dir.z;

    // add a correction factor to the number of absorption lengths
    // abs_lens_left before the photon is absorbed. This factor will be
    // taken out after this propagation step. Usually the factor is 1
    // and thus has no effect, but it is used in a direction-dependent
    // way for our model of ice anisotropy.
    const float abs_len_correction_factor = getDirectionalAbsLenCorrFactor(ph.dir);
    ph.absLength *= abs_len_correction_factor;

    // the "next" medium boundary (either top or bottom, depending on
    // step direction)
    float mediumBoundary = (photon_dz < 0.0f)
                                ? (mediumLayerBoundary(currentPhotonLayer))
                                : (mediumLayerBoundary(currentPhotonLayer) + (float)MEDIUM_LAYER_THICKNESS);

     // track this thing to the next scattering point
    float scaStepLeft = -logf(rng.randUniformFloatOC());

    float currentScaLen = getScatteringLength(currentPhotonLayer, ph.wlen, scatteringLength);
    float currentAbsLen = getAbsorptionLength(currentPhotonLayer, ph.wlen, absorptionDust, absorptionTauDelta);

    float ais = (photon_dz * scaStepLeft - ((mediumBoundary - effective_z)) / currentScaLen) *
                (1.0f / (float)MEDIUM_LAYER_THICKNESS);
    float aia = (photon_dz * ph.absLength - ((mediumBoundary - effective_z)) / currentAbsLen) *
                (1.0f / (float)MEDIUM_LAYER_THICKNESS);

    
    // // propagate through layers
    int j = currentPhotonLayer;
    if (photon_dz < 0) {
        for (; (j > 0) && (ais < 0.0f) && (aia < 0.0f);
                mediumBoundary -= (float)MEDIUM_LAYER_THICKNESS,
                currentScaLen = getScatteringLength(j, ph.wlen, scatteringLength),
                currentAbsLen = getAbsorptionLength(j, ph.wlen, absorptionDust, absorptionTauDelta), 
                ais += 1.f / (currentScaLen),
                aia += 1.f / (currentAbsLen))
            --j;
    } else {
        for (; (j < MEDIUM_LAYERS - 1) && (ais > 0.0f) && (aia > 0.0f);
                mediumBoundary += (float)MEDIUM_LAYER_THICKNESS,
                currentScaLen = getScatteringLength(j, ph.wlen, scatteringLength),
                currentAbsLen = getAbsorptionLength(j, ph.wlen, absorptionDust, absorptionTauDelta), 
                ais -= 1.f / (currentScaLen),
                aia -= 1.f / (currentAbsLen))
            ++j;
    }

    float distanceToAbsorption;
    if ((currentPhotonLayer == j) || ((fabs(photon_dz)) < EPSILON)) {
        distancePropagated = scaStepLeft * currentScaLen;
        distanceToAbsorption = ph.absLength * currentAbsLen;
    } else {
        const float recip_photon_dz = 1.f / (photon_dz);
        distancePropagated =
            (ais * ((float)MEDIUM_LAYER_THICKNESS) * currentScaLen + mediumBoundary - effective_z) *
            recip_photon_dz;
        distanceToAbsorption =
            (aia * ((float)MEDIUM_LAYER_THICKNESS) * currentAbsLen + mediumBoundary - effective_z) *
            recip_photon_dz;
    }

    // get overburden for distance i.e. check if photon is absorbed
    if (distanceToAbsorption < distancePropagated) {
        distancePropagated = distanceToAbsorption;
        ph.absLength = 0.0f;
        return true;
    } else {
        ph.absLength = (distanceToAbsorption - distancePropagated) / currentAbsLen;
        
        // hoist the correction factor back out of the absorption length
        ph.absLength = ph.absLength / abs_len_correction_factor;
        return false;
    }
}

__device__ __forceinline__  void updatePhotonTrack(I3CLPhoton& ph, float distancePropagated)
{
        ph.pos += ph.dir * distancePropagated;
        ph.time += ph.invGroupvel * distancePropagated;
}

__device__ __forceinline__  void scatterPhoton(I3CLPhoton& ph, RngType& rng)
{
     // optional direction transformation (for ice anisotropy)
    detail::transformDirectionPreScatter(ph.dir);

    // choose a scattering angle
    const float cosScatAngle = detail::makeScatteringCosAngle(rng.randUniformFloatCO());
    const float sinScatAngle = sqrt(1.0f - sqr(cosScatAngle));

    // change the current direction by that angle
    detail::scatterDirectionByAngle(cosScatAngle, sinScatAngle, ph.dir, rng.randUniformFloatCO());

    // optional direction transformation (for ice anisotropy)
    detail::transformDirectionPostScatter(ph.dir);
}

#endif