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
#include <cooperative_groups.h>
namespace cg = cooperative_groups;

#include "settings.cuh"
#include "dataStructCuda.cuh"
#include "utils.cuh"
#include "rng.cuh"
#include "domGeoData.cuh"
#include "wlenBiasSource.cuh"
#include "zOffsetHandling.cuh"
#include "wlenGeneration.cuh"
#include "scatteringAndAbsorbtionData.cuh"
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
__device__ __forceinline__ I3CLPhoton createPhoton(const I3CLSimStepCuda& step, const float3& stepDir, const float* wlenLut, RngType& rng);

/**
 * @brief  propgates a single photon once
 * @param ph the photon to propagate
 * @param distanceTraveled the distance the photon was propagated during this iteration
 * @param rng the rng to use for generating this photon, 1 rng values is computed
 * @param scatteringLength scattering length look up table, which can be in global or shared memory
 * @param absorptionDust dust absorption look up table, which can be in global or shared memory
 * @param absorptionTauDelta absorption Tau-Delta look up table, which can be in global or shared memory
 * @param zOffsetLut lut containing zOffset values 
 * @return true if the photon was absorbed
 */
__device__ __forceinline__ bool propPhoton(I3CLPhoton& ph, float& distanceTraveled, RngType& rng, const float* scatteringLength, const float* absorptionDust, const float* absorptionDeltaTau, const float* zOffsetLut);

/**
 * @brief moves a photon along its track by the propagated distance
 * @param ph the photon to move
 * @param distanceTraveled the distance the photon was propagated this iteration
 */
__device__ __forceinline__  void updatePhotonTrack(I3CLPhoton& ph, float distanceTraveled);

/**
 * @brief scatters a photon
 * @param ph the photon to scatter
 * @param rng the random number generator 
 */
__device__ __forceinline__  void scatterPhoton(I3CLPhoton& ph, RngType& rng);

/**
 * @brief perform collision check using the original clsim method,
 *          photons are stored in outputPhotons if a collision is detected
 * 
 */
__device__ __forceinline__ bool checkForCollisionOld(const I3CLPhoton& photon, const I3CLSimStepCuda &step, 
                                                  float &traveledDistance,
                                                  uint32_t *hitIndex, uint32_t maxHitIndex,
                                                  I3CLSimPhotonCuda *outputPhotons,  
                                                  const unsigned short *geoLayerToOMNumIndexPerStringSetLocal,
                                                  const float* getWavelengthBias_dataShared);

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
        dir = float3{(1.0108098679e+00f * dir.x) + (-5.7549125949e-02f * dir.y) + (0.f * dir.z),
                    (-5.7549125949e-02f * dir.x) + (1.0311047952e+00f * dir.y) + (0.f * dir.z),
                    (0.f * dir.x) + (0.f * dir.y) + (9.6252041756e-01f * dir.z)};
        dir = normalize(dir);
    }

    // compute the scattering cosine by selecting between variant 1 and 2
    __device__ __forceinline__ float makeScatteringCosAngle(float randomNumberCO)
    {
        if (randomNumberCO < 3.5000000000e-01f) {
            // const float g = 9.0000000000e-01f;
            // const float beta = (1.f-g)/(1.f+g);
            randomNumberCO = randomNumberCO / 3.5000000000e-01f;
            const float beta = 5.2631578947e-02f;
            return clamp(2.f * powf((randomNumberCO), beta) - 1.f, -1.f, 1.f);
        } else {
            randomNumberCO = (1.f - randomNumberCO) / 6.5000000000e-01f;
            const float g = 9.0000000000e-01f;
            const float g2 = 8.1000000000e-01f;

            // a random number [-1;+1]
            const float s = 2.f * (randomNumberCO)-1.f;

            const float ii = ((1.f - g2) / (1.f + g * s));
            return clamp((1.f + g2 - ii * ii) / (2.f * g), -1.f, 1.f);
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

    // retruns the ice layer number for the z position
    __device__ __forceinline__ int findLayerForGivenZPos(float posZ)
    {
        return int((posZ - (float)MEDIUM_LAYER_BOTTOM_POS) / (float)MEDIUM_LAYER_THICKNESS);
    }

    // returns the z value of the next ice layer boundary
    __device__ __forceinline__ float mediumLayerBoundary(int layer)
    {
        return (float(layer) * ((float)MEDIUM_LAYER_THICKNESS)) + (float)MEDIUM_LAYER_BOTTOM_POS;
    }

    // collision detection helper functions
    // --------------------

    // checks collision between a particular photon and a string of doms
    __device__ __forceinline__ void checkForCollision_OnStringOld(const unsigned short stringNum,
                                                           const float photonDirLenXYSqr,
                                                           const float3 &photonPos,
                                                           const float3 &photonDir,
                                                           float &thisStepLength, bool &hitRecorded,
                                                           unsigned short &hitOnString, unsigned short &hitOnDom,
                                                           const unsigned short *geoLayerToOMNumIndexPerStringSetLocal,
                                                           const unsigned char* geoStringInSet, 
                                                           const unsigned short* geoLayerNumLocal, 
                                                           const float* geoLayerStartZLocal, 
                                                           const float* geoLayerHeightLocal,
                                                           const float* geoStringPosXPointer, 
                                                           const float* geoStringPosYPointer, 
                                                           const float* geoStringMinZPointer, 
                                                           const float* geoStringMaxZPointer)
    {
        // find the string set for this string
        unsigned char stringSet = geoStringInSet[stringNum];

        {  // check intersection with string cylinder
            // only use test if uhat lateral component is bigger than about 0.1 (NEED to
            // check bigger than zero)
            const float smin = sqr(((photonPos.x - float(geoStringPosXPointer[stringNum])) * photonDir.y -
                                    (photonPos.y - float(geoStringPosYPointer[stringNum])) * photonDir.x)) /
                            photonDirLenXYSqr;
            // if (smin > sqr( float(geoStringRadius[stringNum]))) return;  //
            // NOTE: smin == distance squared

            if (smin > sqr(float(GEO_STRING_MAX_RADIUS))) return;  // NOTE: smin == distance squared
        }

        {  // check if photon is above or below the string (geoStringMaxZ and
            // geoStringMinZ do not include the OM radius!)
            if ((photonDir.z > 0.0f) && (photonPos.z > geoStringMaxZPointer[stringNum] + OM_RADIUS)) return;
            if ((photonDir.z < 0.0f) && (photonPos.z < geoStringMinZPointer[stringNum] - OM_RADIUS)) return;
        }

        // this photon could potentially be hitting an om
        // -> check them all

        int lowLayerZ = int((photonPos.z - geoLayerStartZLocal[stringSet]) / geoLayerHeightLocal[stringSet]);
        int highLayerZ = int((photonPos.z + photonDir.z * (thisStepLength)-geoLayerStartZLocal[stringSet]) /
                            geoLayerHeightLocal[stringSet]);
        if (highLayerZ < lowLayerZ) {
            int tmp = lowLayerZ;
            lowLayerZ = highLayerZ;
            highLayerZ = tmp;
        }
        lowLayerZ = min(max(lowLayerZ, 0), geoLayerNumLocal[stringSet] - 1);
        highLayerZ = min(max(highLayerZ, 0), geoLayerNumLocal[stringSet] - 1);

        const unsigned short *geoLayerToOMNumIndex =
            geoLayerToOMNumIndexPerStringSetLocal + (uint32_t(stringSet) * GEO_LAYER_STRINGSET_MAX_NUM_LAYERS) + lowLayerZ;

        for (int layer_z = lowLayerZ; layer_z <= highLayerZ; ++layer_z, ++geoLayerToOMNumIndex) {
            const unsigned short domNum = *geoLayerToOMNumIndex;
            if (domNum == 0xFFFF) continue;  // empty layer for this string

            float domPosX, domPosY, domPosZ;
            geometryGetDomPosition(stringNum, domNum, domPosX, domPosY, domPosZ);

            float urdot, discr;
            {
                const float4 drvec =
                    float4{domPosX - photonPos.x, domPosY - photonPos.y, domPosZ - photonPos.z, 0.0f};
                const float dr2 = dot(drvec, drvec);

                urdot = dot(drvec, make_float4(photonDir,0));              // this assumes drvec.w==0
                discr = sqr(urdot) - dr2 + OM_RADIUS * OM_RADIUS;  // (discr)^2
            }

            if (discr < 0.0f) continue;  // no intersection with this DOM
            discr = sqrtf(discr);

            // by construction: smin1 < smin2
            {
                // distance from current point along the track to second intersection
                const float smin2 = urdot + discr;

                if (smin2 < 0.0f) continue;  // implies smin1 < 0, so no intersection
            }

            // distance from current point along the track to first intersection
            const float smin1 = urdot - discr;

            // smin2 > 0 && smin1 < 0 means that there *is* an intersection, but we are
            // starting inside the DOM. This allows photons starting inside a DOM to
            // leave (necessary for flashers):
            if (smin1 < 0.0f) continue;

                // if we get here, there *is* an intersection with the DOM (there are two
                // actually, one for the ray enetering the DOM and one when it leaves
                // again). We are interested in the one where enters the ray enters the
                // DOM.

                // check if distance to intersection <= thisStepLength; if not then no
                // detection
            if (smin1 < thisStepLength)
            {
                // record a hit (for later, the actual recording is done
                // in checkForCollision().)
                thisStepLength = smin1;  // limit step length
                hitOnString = stringNum;
                hitOnDom = domNum;
                hitRecorded = true;
                // continue searching, maybe we hit a closer OM..
                // (in that case, no hit will be saved for this one)
            }
        }  // end for loop layer_z
    }

    // calls checkForCollision_OnStringOld() for every string that a photon might collide with based on the x/y grid
    __device__ __forceinline__ void checkForCollision_InCellOld(
        const float photonDirLenXYSqr, const float3 &photonPos, const float3 &photonDir,
        float &thisStepLength, bool &hitRecorded, unsigned short &hitOnString, unsigned short &hitOnDom,
        const unsigned short *geoLayerToOMNumIndexPerStringSetLocal, const unsigned short *this_geoCellIndex,
        const float this_geoCellStartX, const float this_geoCellStartY, const float this_geoCellWidthX,
        const float this_geoCellWidthY, const int this_geoCellNumX, const int this_geoCellNumY,
        const unsigned char* geoStringInSet, 
        const unsigned short* geoLayerNumLocal, 
        const float* geoLayerStartZLocal, 
        const float* geoLayerHeightLocal,
        const float* geoStringPosXPointer, 
        const float* geoStringPosYPointer,
        const float* geoStringMinZPointer, 
        const float* geoStringMaxZPointer)
    {
        int lowCellX = int((photonPos.x - this_geoCellStartX) / this_geoCellWidthX);
        int lowCellY = int((photonPos.y - this_geoCellStartY) / this_geoCellWidthY);

        int highCellX =
            int((photonPos.x + photonDir.x * (thisStepLength)-this_geoCellStartX) / this_geoCellWidthX);
        int highCellY =
            int((photonPos.y + photonDir.y * (thisStepLength)-this_geoCellStartY) / this_geoCellWidthY);

        if (highCellX < lowCellX) {
            int tmp = lowCellX;
            lowCellX = highCellX;
            highCellX = tmp;
        }
        if (highCellY < lowCellY) {
            int tmp = lowCellY;
            lowCellY = highCellY;
            highCellY = tmp;
        }

        lowCellX = min(max(lowCellX, 0), this_geoCellNumX - 1);
        lowCellY = min(max(lowCellY, 0), this_geoCellNumY - 1);
        highCellX = min(max(highCellX, 0), this_geoCellNumX - 1);
        highCellY = min(max(highCellY, 0), this_geoCellNumY - 1);

        for (int cell_y = lowCellY; cell_y <= highCellY; ++cell_y) {
            for (int cell_x = lowCellX; cell_x <= highCellX; ++cell_x) {
                const unsigned short stringNum = this_geoCellIndex[cell_y * this_geoCellNumX + cell_x];
                if (stringNum == 0xFFFF) continue;  // empty cell
                checkForCollision_OnStringOld(stringNum, photonDirLenXYSqr, photonPos, photonDir, 
                                        thisStepLength, hitRecorded, hitOnString, hitOnDom,
                                        geoLayerToOMNumIndexPerStringSetLocal, 
                                        geoStringInSet, geoLayerNumLocal, geoLayerStartZLocal, geoLayerHeightLocal,
                                        geoStringPosXPointer, geoStringPosYPointer, geoStringMinZPointer, geoStringMaxZPointer);
            }
        }
    }


    // calls checkForCollision_InCellOld() for every individual x/y grid that is defined (up to two!)
    __device__ __forceinline__ void checkForCollision_InCellsOld(const float photonDirLenXYSqr, const float3 &photonPos,
                                                          const float3 &photonDir,
                                                          float &thisStepLength, bool &hitRecorded,
                                                          unsigned short &hitOnString, unsigned short &hitOnDom,
                                                          const unsigned short *geoLayerToOMNumIndexPerStringSetLocal,
                                                          const unsigned short* geoCellIndex0, 
                                                          const unsigned short* geoCellIndex1, 
                                                          const unsigned char* geoStringInSet, 
                                                          const unsigned short* geoLayerNumLocal, 
                                                          const float* geoLayerStartZLocal, 
                                                          const float* geoLayerHeightLocal,
                                                          const float* geoStringPosXPointer, 
                                                          const float* geoStringPosYPointer,
                                                          const float* geoStringMinZPointer, 
                                                          const float* geoStringMaxZPointer)
    {
        // using macros and hard-coded names is
        // not really the best thing to do here..
        // replace with a loop sometime.
    #define DO_CHECK(subdetectorNum)                                                                                 \
        checkForCollision_InCellOld(photonDirLenXYSqr, photonPos, photonDir, thisStepLength, hitRecorded, \
                                hitOnString, hitOnDom, geoLayerToOMNumIndexPerStringSetLocal,                       \
                                                                                                                    \
                                geoCellIndex##subdetectorNum, GEO_CELL_START_X_##subdetectorNum,  \
                                GEO_CELL_START_Y_##subdetectorNum, GEO_CELL_WIDTH_X_##subdetectorNum,               \
                                GEO_CELL_WIDTH_Y_##subdetectorNum, GEO_CELL_NUM_X_##subdetectorNum,                 \
                                GEO_CELL_NUM_Y_##subdetectorNum,                                                    \
                                geoStringInSet, geoLayerNumLocal, geoLayerStartZLocal, geoLayerHeightLocal,         \
                                geoStringPosXPointer, geoStringPosYPointer, geoStringMinZPointer, geoStringMaxZPointer);        

        // argh..
    #if GEO_CELL_NUM_SUBDETECTORS > 0
        DO_CHECK(0);
    #endif

    #if GEO_CELL_NUM_SUBDETECTORS > 1
        DO_CHECK(1);
    #endif
    #undef DO_CHECK

        // TODO: add additional defines for future geometry or bettwe, use different collision detection algorithm 
    }

    // calculate spherical coordinates from cartesian coordinates to store hit point
    __device__ __forceinline__ float2 sphDirFromCar(float3 carDir)
    {
        // Calculate Spherical coordinates from Cartesian
        const float r_inv = rsqrtf(carDir.x * carDir.x + carDir.y * carDir.y + carDir.z * carDir.z);

        float theta = 0.f;
        if ((fabs(carDir.z * r_inv)) <= 1.f) {
            theta = acos(carDir.z * r_inv);
        } else {
            if (carDir.z < 0.f) theta = PI;
        }
        if (theta < 0.f) theta += 2.f * PI;

        float phi = atan2(carDir.y, carDir.x);
        if (phi < 0.f) phi += 2.f * PI;

        return float2{theta, phi};
    }

    // stores a hit in outputPhotons as long as there is still room
    __device__ __forceinline__ void saveHit( const I3CLPhoton& photon, float traveledDistance,
                                                unsigned short hitOnString, unsigned short hitOnDom,
                                                uint32_t *hitIndex, uint32_t maxHitIndex, I3CLSimPhotonCuda *outputPhotons,
                                                const float* getWavelengthBias_dataShared,
                                                const I3CLSimStepCuda &step)
    {
        uint32_t myIndex = atomicAdd(hitIndex, 1);

        if (myIndex < maxHitIndex) {
            // Emit photon position relative to the hit DOM

            float domPosX, domPosY, domPosZ;
            geometryGetDomPosition(hitOnString, hitOnDom, domPosX, domPosY, domPosZ);

            I3CLSimPhotonCuda outphoton;
            outphoton.posAndTime = float4{photon.pos.x + traveledDistance * photon.dir.x - domPosX,
                                        photon.pos.y + traveledDistance * photon.dir.y - domPosY,
                                        photon.pos.z + traveledDistance * photon.dir.z - domPosZ,
                                        photon.time + traveledDistance * photon.invGroupvel};
            outphoton.dir = sphDirFromCar(photon.dir);
            outphoton.wavelength = photon.wlen;
            outphoton.weight = step.weight / getWavelengthBias(photon.wlen, getWavelengthBias_dataShared);
            outphoton.groupVelocity = 1.f / (photon.invGroupvel);
            outphoton.identifier = step.identifier;
            outphoton.stringID = hitOnString;
            outphoton.omID = hitOnDom;

            outputPhotons[myIndex] = outphoton;
        }
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

__device__ __forceinline__ I3CLPhoton createPhoton(const I3CLSimStepCuda& step, const float3& stepDir, const float* wlenLut, RngType& rng)
{
    // float4 randomNumbers = float4{0.12454854,0.99568,0.4877858,0.24784564};
    float4 randomNumbers = {rng.randUniformFloatCO(), rng.randUniformFloatCO(), rng.randUniformFloatCO(), rng.randUniformFloatOC()};

    I3CLPhoton ph;

    // move along the step direction a random amount
    const float shiftMultiplied = step.dirAndLengthAndBeta.z * randomNumbers.x;
    const float inverseParticleSpeed = 1.f / (C_LIGHT * step.dirAndLengthAndBeta.w);
    ph.pos = float3{step.posAndTime.x + stepDir.x * shiftMultiplied, step.posAndTime.y + stepDir.y * shiftMultiplied,
                    step.posAndTime.z + stepDir.z * shiftMultiplied};
    ph.time = step.posAndTime.w + inverseParticleSpeed * shiftMultiplied;

    // generate a wavelength
    ph.wlen = getWavelenth(randomNumbers.y,wlenLut);

    // calculate phase and group ref index
    const float x = ph.wlen * 1e6;
    const float phaseRefIndex = detail::getPhaseRefIndex(ph.wlen,x);
    const float groupRefIndex = detail::getGroupRefIndex(ph.wlen,phaseRefIndex,x);

    const float cosCherenkov = min( 1.0f, 1.0f / (step.dirAndLengthAndBeta.w * phaseRefIndex));  // cos theta = 1/(beta*n)
    const float sinCherenkov = sqrtf(1.0f - sqr(cosCherenkov));
    
    // determine the photon direction
    // start with the track direction and rotate to cherenkov emission direction
    ph.dir = stepDir;
    detail::scatterDirectionByAngle(cosCherenkov, sinCherenkov, ph.dir, randomNumbers.z);
    
    // calc inverse group velocity and set an initial absorption length
    ph.invGroupvel = groupRefIndex * RECIP_C_LIGHT; // refIndex * (1/c_light) <=> 1 / (c_light / refIndex)
    ph.absLength = -logf(randomNumbers.w);

    return ph;
}

__device__ __forceinline__ bool propPhoton(I3CLPhoton& ph, float& distancePropagated, RngType& rng, const float* scatteringLength, const float* absorptionDust, const float* absorptionTauDelta, const float* zOffsetLut)
{ 
    const float effective_z = ph.pos.z - getZOffset(ph.pos, zOffsetLut); //ref::getTiltZShift(make_float4(ph.pos,0)); 
    const int currentPhotonLayer = min(max( detail::findLayerForGivenZPos(effective_z), 0), MEDIUM_LAYERS - 1);
    const float photon_dz = ph.dir.z;

    // add a correction factor to the number of absorption lengths
    // abs_lens_left before the photon is absorbed. This factor will be
    // taken out after this propagation step. Usually the factor is 1
    // and thus has no effect, but it is used in a direction-dependent
    // way for our model of ice anisotropy.
    const float abs_len_correction_factor = detail::getDirectionalAbsLenCorrFactor(ph.dir);
    ph.absLength *= abs_len_correction_factor;

    // the "next" medium boundary (either top or bottom, depending on
    // step direction)
    float mediumBoundary = (photon_dz < 0.0f)
                                ? (detail::mediumLayerBoundary(currentPhotonLayer))
                                : (detail::mediumLayerBoundary(currentPhotonLayer) + (float)MEDIUM_LAYER_THICKNESS);

    // track this thing to the next scattering point
    float scaStepLeft;
    #ifdef BLOCK_RANDOM_NUMBERS_PROPAGATION
        cg::coalesced_group active = cg::coalesced_threads();
        if(active.thread_rank()==0)
            scaStepLeft = -logf(rng.randUniformFloatOC());
        scaStepLeft = active.shfl(scaStepLeft, 0);
    #else
        scaStepLeft = -logf(rng.randUniformFloatOC());
    #endif

    float currentScaLen = detail::getScatteringLength(currentPhotonLayer, ph.wlen, scatteringLength);
    float currentAbsLen = detail::getAbsorptionLength(currentPhotonLayer, ph.wlen, absorptionDust, absorptionTauDelta);

    float ais = (photon_dz * scaStepLeft - ((mediumBoundary - effective_z)) / currentScaLen) *
                (1.0f / (float)MEDIUM_LAYER_THICKNESS);
    float aia = (photon_dz * ph.absLength - ((mediumBoundary - effective_z)) / currentAbsLen) *
                (1.0f / (float)MEDIUM_LAYER_THICKNESS);

    
    float dir = copysign(1.0f,photon_dz);
    int j = currentPhotonLayer;
    for (; (j > 0) && (j < MEDIUM_LAYERS - 1) && (ais*dir > 0.0f) && (aia*dir > 0.0f);
            mediumBoundary += dir* (float)MEDIUM_LAYER_THICKNESS,
            currentScaLen = detail::getScatteringLength(j, ph.wlen, scatteringLength),
            currentAbsLen = detail::getAbsorptionLength(j, ph.wlen, absorptionDust, absorptionTauDelta), 
            ais -= dir / (currentScaLen),
            aia -= dir / (currentAbsLen))
        --j;

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

    float r1;
    float r2;
    #ifdef BLOCK_RANDOM_NUMBERS_PROPAGATION
        cg::coalesced_group active = cg::coalesced_threads();
        if(active.thread_rank()==0)
        {
            r1 = rng.randUniformFloatCO();
            r2 = rng.randUniformFloatCO();
        }
        r1 = active.shfl(r1,0);
        r2 = active.shfl(r2,0);
    #else
        r1 = rng.randUniformFloatCO();
        r2 = rng.randUniformFloatCO();
    #endif

    // choose a scattering angle
    const float cosScatAngle = detail::makeScatteringCosAngle(r1);
    const float sinScatAngle = sqrt(1.0f - sqr(cosScatAngle));

    // change the current direction by that angle
    detail::scatterDirectionByAngle(cosScatAngle, sinScatAngle, ph.dir, r2);

    // optional direction transformation (for ice anisotropy)
    detail::transformDirectionPostScatter(ph.dir);
}

__device__ __forceinline__ bool checkForCollisionOld(const I3CLPhoton& photon, const I3CLSimStepCuda &step, 
                                                  float &traveledDistance,
                                                  uint32_t *hitIndex, uint32_t maxHitIndex,
                                                  I3CLSimPhotonCuda *outputPhotons,  
                                                  const unsigned short *geoLayerToOMNumIndexPerStringSetLocal,
                                                  const float* getWavelengthBias_dataShared,
                                                  const unsigned short* geoCellIndex0, 
                                                  const unsigned short* geoCellIndex1, 
                                                  const unsigned char* geoStringInSet, 
                                                  const unsigned short* geoLayerNumLocal, 
                                                  const float* geoLayerStartZLocal, 
                                                  const float* geoLayerHeightLocal,
                                                  const float* geoStringPosXPointer, 
                                                  const float* geoStringPosYPointer,  
                                                  const float* geoStringMinZPointer, 
                                                  const float* geoStringMaxZPointer)
{
    const float photonDirLenXYSqr = sqr(photon.dir.x) + sqr(photon.dir.y); 
    if (photonDirLenXYSqr <= 0.0f) return false; // is this really correct? what if we are inbetween doms and move straight down 

    bool hitRecorded = false;
    unsigned short hitOnString;
    unsigned short hitOnDom;

    detail::checkForCollision_InCellsOld(photonDirLenXYSqr, photon.pos, photon.dir,
                              traveledDistance, hitRecorded, hitOnString, hitOnDom,
                              geoLayerToOMNumIndexPerStringSetLocal, geoCellIndex0, geoCellIndex1, geoStringInSet, 
                              geoLayerNumLocal, geoLayerStartZLocal, geoLayerHeightLocal,
                              geoStringPosXPointer, geoStringPosYPointer, geoStringMinZPointer, geoStringMaxZPointer);

    // record the photons, i.e. store the hit in the outputPhotons array
    if (hitRecorded) {
        detail::saveHit(photon, traveledDistance, hitOnString, hitOnDom, hitIndex, maxHitIndex, outputPhotons, getWavelengthBias_dataShared, step);
    }
    return hitRecorded;
}


#endif