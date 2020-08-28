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

#ifndef PROPAGATIONKERNELFUNCTIUONS_CUH
#define PROPAGATIONKERNELFUNCTIUONS_CUH

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <math_constants.h>

#include <mwcrngKernelSource.cuh>
#include <CLnoneCUDA.cuh>
#include <prependSource.cuh>
#include <wlenGeneratorSource.cuh>
#include <dataStructCuda.cuh>
#include <wlenBiasSource.cuh>
#include <mediumPropertiesSource.cuh>
#include <geometrySource.cuh>

#include <chrono>
#include <ctime>

// for debugging:
#define PRINTLROOT                                                              \
    if (blockIdx.x * blockDim.x + threadIdx.x == 0) {                           \
        printf("thread 0 - in line %d and function %s \n", __LINE__, __func__); \
    }
//#define PRINTL        printf("thread %d - in line %d and function %s \n", blockIdx.x * blockDim.x + threadIdx.x,
//__LINE__, __func__);

///////////////// forward declarations

__device__ __forceinline__ int findLayerForGivenZPos(float posZ);

__device__ __forceinline__ float mediumLayerBoundary(int layer);

__device__ __forceinline__ void scatterDirectionByAngle(float cosa, float sina, float4 &direction, float randomNumber);

__device__ __forceinline__ void createPhotonFromTrack(const I3CLSimStepCuda &step, const float4 &stepDir, RNG_ARGS,
                                                      float4 &photonPosAndTime, float4 &photonDirAndWlen,
                                                      float* _generateWavelength_0distYValuesShared, float* _generateWavelength_0distYCumulativeValuesShared);

__device__ __forceinline__ float2 sphDirFromCar(float4 carDir);

__device__ __forceinline__ void saveHit(const float4 &photonPosAndTime, const float4 &photonDirAndWlen, const float* getWavelengthBias_dataShared,
                                        const float thisStepLength, float inv_groupvel, float photonTotalPathLength,
                                        uint32_t photonNumScatters, float distanceTraveledInAbsorptionLengths,
                                        const float4 &photonStartPosAndTime, const float4 &photonStartDirAndWlen,
                                        const I3CLSimStepCuda &step, unsigned short hitOnString,
                                        unsigned short hitOnDom, uint32_t *hitIndex, uint32_t maxHitIndex,
                                        I3CLSimPhotonCuda *outputPhotons
#ifdef SAVE_PHOTON_HISTORY
                                        ,
                                        float4 *photonHistory, float4 *currentPhotonHistory
#endif
);

///////////////////////// some constants

constexpr float speedOfLight = 0.299792458f;       // [m/ns]F
constexpr float recip_speedOfLight = 3.33564095f;  // [ns/m]
constexpr float PI = 3.14159265359f;

///////////////////////////

__device__ __forceinline__ float my_fabs(const float a);
__device__ __forceinline__ float sqr(const float a);

__device__ __forceinline__ void checkForCollision_OnString(const unsigned short stringNum,
                                                           const float photonDirLenXYSqr,
                                                           const float4 &photonPosAndTime,
                                                           const float4 &photonDirAndWlen,
                                                           const float* getWavelengthBias_dataShared,
#ifdef STOP_PHOTONS_ON_DETECTION
                                                           float &thisStepLength, bool &hitRecorded,
                                                           unsigned short &hitOnString, unsigned short &hitOnDom,
#else
                                                           float thisStepLength, float inv_groupvel,
                                                           float photonTotalPathLength, uint32_t photonNumScatters,
                                                           float distanceTraveledInAbsorptionLengths,
                                                           const float4 &photonStartPosAndTime,
                                                           const float4 &photonStartDirAndWlen,
                                                           const I3CLSimStepCuda &step, uint32_t *hitIndex,
                                                           uint32_t maxHitIndex, I3CLSimPhotonCuda *outputPhotons,
#ifdef SAVE_PHOTON_HISTORY
                                                           float4 *photonHistory, float4 *currentPhotonHistory,
#endif
#endif
                                                           const unsigned short *geoLayerToOMNumIndexPerStringSetLocal);

__device__ __forceinline__ void checkForCollision_InCell(const float photonDirLenXYSqr, const float4 &photonPosAndTime,
                                                         const float4 &photonDirAndWlen, const float* getWavelengthBias_dataShared,
#ifdef STOP_PHOTONS_ON_DETECTION
                                                         float &thisStepLength, bool &hitRecorded,
                                                         unsigned short &hitOnString, unsigned short &hitOnDom,
#else
                                                         float thisStepLength, float inv_groupvel,
                                                         float photonTotalPathLength, uint32_t photonNumScatters,
                                                         float distanceTraveledInAbsorptionLengths,
                                                         const float4 &photonStartPosAndTime,
                                                         const float4 &photonStartDirAndWlen,
                                                         const I3CLSimStepCuda &step, uint32_t *hitIndex,
                                                         uint32_t maxHitIndex, I3CLSimPhotonCuda *outputPhotons, 
#ifdef SAVE_PHOTON_HISTORY
                                                         float4 *photonHistory, float4 *currentPhotonHistory,
#endif
#endif
                                                         const unsigned short *geoLayerToOMNumIndexPerStringSetLocal,

                                                         unsigned short *this_geoCellIndex,
                                                         const float this_geoCellStartX, const float this_geoCellStartY,
                                                         const float this_geoCellWidthX, const float this_geoCellWidthY,
                                                         const int this_geoCellNumX, const int this_geoCellNumY);

__device__ __forceinline__ void checkForCollision_InCells(const float photonDirLenXYSqr, const float4 &photonPosAndTime,
                                                          const float4 &photonDirAndWlen, const float* getWavelengthBias_dataShared,
#ifdef STOP_PHOTONS_ON_DETECTION
                                                          float &thisStepLength, bool &hitRecorded,
                                                          unsigned short &hitOnString, unsigned short &hitOnDom,
#else
                                                          float thisStepLength, float inv_groupvel,
                                                          float photonTotalPathLength, uint32_t photonNumScatters,
                                                          float distanceTraveledInAbsorptionLengths,
                                                          const float4 &photonStartPosAndTime,
                                                          const float4 &photonStartDirAndWlen,
                                                          const I3CLSimStepCuda &step, uint32_t *hitIndex,
                                                          uint32_t maxHitIndex, I3CLSimPhotonCuda *outputPhotons, 
#ifdef SAVE_PHOTON_HISTORY
                                                          float4 *photonHistory, float4 *currentPhotonHistory,
#endif
#endif
                                                          const unsigned short *geoLayerToOMNumIndexPerStringSetLocal);

__device__ __forceinline__ bool checkForCollision(const I3CLPhoton& photon, const I3CLInitialPhoton& photonInitials,
                                                  const I3CLSimStepCuda &step, float &thisStepLength,
                                                  uint32_t *hitIndex, uint32_t maxHitIndex,
                                                  I3CLSimPhotonCuda *outputPhotons,  
                                                  const unsigned short *geoLayerToOMNumIndexPerStringSetLocal, 
                                                  const float* getWavelengthBias_dataShared, int* collisionCheckIndex, CollisionCheck* collisionChecks);

__device__ __forceinline__ void checkForCollision_OnString(const unsigned short stringNum,
                                                           const float photonDirLenXYSqr,
                                                           const float4 &photonPosAndTime,
                                                           const float4 &photonDirAndWlen, const float* getWavelengthBias_dataShared,
#ifdef STOP_PHOTONS_ON_DETECTION
                                                           float &thisStepLength, bool &hitRecorded,
                                                           unsigned short &hitOnString, unsigned short &hitOnDom,
#else
                                                           float thisStepLength, float inv_groupvel,
                                                           float photonTotalPathLength, uint32_t photonNumScatters,
                                                           float distanceTraveledInAbsorptionLengths,
                                                           const float4 &photonStartPosAndTime,
                                                           const float4 &photonStartDirAndWlen,
                                                           const I3CLSimStepCuda &step, uint32_t *hitIndex,
                                                           uint32_t maxHitIndex, I3CLSimPhotonCuda *outputPhotons,  
#ifdef SAVE_PHOTON_HISTORY
                                                           float4 *photonHistory, float4 *currentPhotonHistory,
#endif
#endif
                                                           const unsigned short *geoLayerToOMNumIndexPerStringSetLocal)
{
    // find the string set for this string
    unsigned char stringSet = geoStringInStringSet[stringNum];

    {  // check intersection with string cylinder
        // only use test if uhat lateral component is bigger than about 0.1 (NEED to
        // check bigger than zero)
        const float smin = sqr(((photonPosAndTime.x - float(geoStringPosX[stringNum])) * photonDirAndWlen.y -
                                (photonPosAndTime.y - float(geoStringPosY[stringNum])) * photonDirAndWlen.x)) /
                           photonDirLenXYSqr;
        // if (smin > sqr( float(geoStringRadius[stringNum]))) return;  //
        // NOTE: smin == distance squared

        if (smin > sqr(float(GEO_STRING_MAX_RADIUS))) return;  // NOTE: smin == distance squared
    }

    {  // check if photon is above or below the string (geoStringMaxZ and
        // geoStringMinZ do not include the OM radius!)
        if ((photonDirAndWlen.z > ZERO) && (photonPosAndTime.z > geoStringMaxZ[stringNum] + OM_RADIUS)) return;
        if ((photonDirAndWlen.z < ZERO) && (photonPosAndTime.z < geoStringMinZ[stringNum] - OM_RADIUS)) return;
    }

    // this photon could potentially be hitting an om
    // -> check them all

    int lowLayerZ = int((photonPosAndTime.z - geoLayerStartZ[stringSet]) / geoLayerHeight[stringSet]);
#ifdef STOP_PHOTONS_ON_DETECTION
    int highLayerZ = int((photonPosAndTime.z + photonDirAndWlen.z * (thisStepLength)-geoLayerStartZ[stringSet]) /
                         geoLayerHeight[stringSet]);
#else
    int highLayerZ = int((photonPosAndTime.z + photonDirAndWlen.z * thisStepLength - geoLayerStartZ[stringSet]) /
                         geoLayerHeight[stringSet]);
#endif
    if (highLayerZ < lowLayerZ) {
        int tmp = lowLayerZ;
        lowLayerZ = highLayerZ;
        highLayerZ = tmp;
    }
    lowLayerZ = min(max(lowLayerZ, 0), geoLayerNum[stringSet] - 1);
    highLayerZ = min(max(highLayerZ, 0), geoLayerNum[stringSet] - 1);

#ifndef STOP_PHOTONS_ON_DETECTION
// the number of 64bit integers needed to store bits for all doms
#define numComponents ((GEO_MAX_DOM_INDEX + 64 - 1) / 64)
    uint64_t dom_bitmask[numComponents];
    for (uint32_t i = 0; i < numComponents; ++i) dom_bitmask[i] = 0;
#undef numComponents
#endif

    //__device__ const unsigned short
    //*geoLayerToOMNumIndex=geoLayerToOMNumIndexPerStringSet + (
    // uint(stringSet)*GEO_LAYER_STRINGSET_MAX_NUM_LAYERS) + lowLayerZ;
    //  __local const unsigned short *geoLayerToOMNumIndex=geoLayerToOMNumIndexPerStringSetLocal +
    //  (convert_uint(stringSet)*GEO_LAYER_STRINGSET_MAX_NUM_LAYERS) + lowLayerZ;
    const unsigned short *geoLayerToOMNumIndex =
        geoLayerToOMNumIndexPerStringSetLocal + (uint32_t(stringSet) * GEO_LAYER_STRINGSET_MAX_NUM_LAYERS) + lowLayerZ;

    for (int layer_z = lowLayerZ; layer_z <= highLayerZ; ++layer_z, ++geoLayerToOMNumIndex) {
        const unsigned short domNum = *geoLayerToOMNumIndex;
        if (domNum == 0xFFFF) continue;  // empty layer for this string

#ifndef STOP_PHOTONS_ON_DETECTION
        // prevent strings from being checked twice
        if (dom_bitmask[stringNum / 64] & (1 << uint64_t(domNum % 64))) continue;  // already check this string
        dom_bitmask[stringNum / 64] |= (1 << uint64_t(domNum % 64));               // mark this string as checked
#endif

#ifndef CABLE_RADIUS
        float domPosX, domPosY, domPosZ;
        geometryGetDomPosition(stringNum, domNum, domPosX, domPosY, domPosZ);
#else
        float domPosX, domPosY, domPosZ, cableOrientation;
        geometryGetDomPosition(stringNum, domNum, domPosX, domPosY, domPosZ, &cableOrientation);
#endif

        float urdot, discr;
        {
            const float4 drvec =
                float4{domPosX - photonPosAndTime.x, domPosY - photonPosAndTime.y, domPosZ - photonPosAndTime.z, ZERO};
            const float dr2 = dot(drvec, drvec);

            urdot = dot(drvec, photonDirAndWlen);              // this assumes drvec.w==0
            discr = sqr(urdot) - dr2 + OM_RADIUS * OM_RADIUS;  // (discr)^2
        }

#ifdef CABLE_RADIUS
        // Check intersection with cable
        float discr_cable;
        {
            // check intersection with infinite cylinder
            const float4 drvec =
                float4{domPosX + (OM_RADIUS + CABLE_RADIUS) * cos(cableOrientation) - photonPosAndTime.x,
                       domPosY + (OM_RADIUS + CABLE_RADIUS) * sin(cableOrientation) - photonPosAndTime.y, ZERO, ZERO};
            const float dr2 = dot(drvec, drvec);

            const float h_norm = hypot(photonDirAndWlen.x, photonDirAndWlen.y);
            const float urdot =
                h_norm > ZERO ? dot(drvec, photonDirAndWlen / h_norm) : ZERO;  // this assumes drvec.w==0
            discr_cable = sqr(urdot) - dr2 + CABLE_RADIUS * CABLE_RADIUS;      // (discr)^2
        }

        // no intersection, or blocked by cable
        if (discr < ZERO || discr_cable >= ZERO) continue;
#else
        if (discr < ZERO) continue;  // no intersection with this DOM
#endif

#ifdef PANCAKE_FACTOR
        discr = sqrtf(discr) / PANCAKE_FACTOR;
#else
        discr = sqrtf(discr);
#endif

        // by construction: smin1 < smin2

        {
            // distance from current point along the track to second intersection
            const float smin2 = urdot + discr;

            if (smin2 < ZERO) continue;  // implies smin1 < 0, so no intersection
        }

        // distance from current point along the track to first intersection
        const float smin1 = urdot - discr;

        // smin2 > 0 && smin1 < 0 means that there *is* an intersection, but we are
        // starting inside the DOM. This allows photons starting inside a DOM to
        // leave (necessary for flashers):
        if (smin1 < ZERO) continue;

            // if we get here, there *is* an intersection with the DOM (there are two
            // actually, one for the ray enetering the DOM and one when it leaves
            // again). We are interested in the one where enters the ray enters the
            // DOM.

            // check if distance to intersection <= thisStepLength; if not then no
            // detection
#ifdef STOP_PHOTONS_ON_DETECTION
        if (smin1 < thisStepLength)
#else
        if (smin1 < thisStepLength)
#endif
        {
#ifdef STOP_PHOTONS_ON_DETECTION
            // record a hit (for later, the actual recording is done
            // in checkForCollision().)
            thisStepLength = smin1;  // limit step length
            hitOnString = stringNum;
            hitOnDom = domNum;
            hitRecorded = true;
            // continue searching, maybe we hit a closer OM..
            // (in that case, no hit will be saved for this one)
#else  // STOP_PHOTONS_ON_DETECTION
       // save the hit right here

            saveHit(photonPosAndTime, photonDirAndWlen, getWavelengthBias_dataShared,
                    smin1,  // this is the limited thisStepLength
                    inv_groupvel, photonTotalPathLength, photonNumScatters, distanceTraveledInAbsorptionLengths,
                    photonStartPosAndTime, photonStartDirAndWlen, step, stringNum, domNum, hitIndex, maxHitIndex,
                    outputPhotons
#ifdef SAVE_PHOTON_HISTORY
                    ,
                    photonHistory, currentPhotonHistory
#endif  // SAVE_PHOTON_HISTORY
            );
#endif  // STOP_PHOTONS_ON_DETECTION
        }
    }  // end for loop layer_z
}

__device__ __forceinline__ void checkForCollision_InCell(
    const float photonDirLenXYSqr, const float4 &photonPosAndTime, const float4 &photonDirAndWlen, const float* getWavelengthBias_dataShared,
#ifdef STOP_PHOTONS_ON_DETECTION
    float &thisStepLength, bool &hitRecorded, unsigned short &hitOnString, unsigned short &hitOnDom,
#else
    float thisStepLength, float inv_groupvel, float photonTotalPathLength, uint32_t photonNumScatters,
    float distanceTraveledInAbsorptionLengths, const float4 &photonStartPosAndTime, const float4 &photonStartDirAndWlen,
    const I3CLSimStepCuda &step, uint32_t *hitIndex, uint32_t maxHitIndex, I3CLSimPhotonCuda *outputPhotons,  
#ifdef SAVE_PHOTON_HISTORY
    float4 *photonHistory, float4 *currentPhotonHistory,
#endif
#endif
    const unsigned short *geoLayerToOMNumIndexPerStringSetLocal, unsigned short *this_geoCellIndex,
    const float this_geoCellStartX, const float this_geoCellStartY, const float this_geoCellWidthX,
    const float this_geoCellWidthY, const int this_geoCellNumX, const int this_geoCellNumY)
{
    int lowCellX = int((photonPosAndTime.x - this_geoCellStartX) / this_geoCellWidthX);
    int lowCellY = int((photonPosAndTime.y - this_geoCellStartY) / this_geoCellWidthY);

#ifdef STOP_PHOTONS_ON_DETECTION
    int highCellX =
        int((photonPosAndTime.x + photonDirAndWlen.x * (thisStepLength)-this_geoCellStartX) / this_geoCellWidthX);
    int highCellY =
        int((photonPosAndTime.y + photonDirAndWlen.y * (thisStepLength)-this_geoCellStartY) / this_geoCellWidthY);
#else
    int highCellX =
        int((photonPosAndTime.x + photonDirAndWlen.x * thisStepLength - this_geoCellStartX) / this_geoCellWidthX);
    int highCellY =
        int((photonPosAndTime.y + photonDirAndWlen.y * thisStepLength - this_geoCellStartY) / this_geoCellWidthY);
#endif

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

#ifndef STOP_PHOTONS_ON_DETECTION
// the number of 64bit integers needed to store bits for all strings
#define numComponents ((NUM_STRINGS + 64 - 1) / 64)
    uint64_t string_bitmask[numComponents];
    for (uint32_t i = 0; i < numComponents; ++i) string_bitmask[i] = 0;
#undef numComponents
#endif

    for (int cell_y = lowCellY; cell_y <= highCellY; ++cell_y) {
        for (int cell_x = lowCellX; cell_x <= highCellX; ++cell_x) {
            const unsigned short stringNum = this_geoCellIndex[cell_y * this_geoCellNumX + cell_x];
            if (stringNum == 0xFFFF) continue;  // empty cell

#ifndef STOP_PHOTONS_ON_DETECTION
            // prevent strings from being checked twice
            if (string_bitmask[stringNum / 64] & (1 << uint64_t(stringNum % 64)))
                continue;                                                       // already check this string
            string_bitmask[stringNum / 64] |= (1 << uint64_t(stringNum % 64));  // mark this string as checked
#endif

            checkForCollision_OnString(stringNum, photonDirLenXYSqr, photonPosAndTime, photonDirAndWlen, getWavelengthBias_dataShared,
#ifdef STOP_PHOTONS_ON_DETECTION
                                       thisStepLength, hitRecorded, hitOnString, hitOnDom,
#else  // STOP_PHOTONS_ON_DETECTION
                                       thisStepLength, inv_groupvel, photonTotalPathLength, photonNumScatters,
                                       distanceTraveledInAbsorptionLengths, photonStartPosAndTime,
                                       photonStartDirAndWlen, step, hitIndex, maxHitIndex, outputPhotons
#ifdef SAVE_PHOTON_HISTORY
                                       photonHistory, currentPhotonHistory,
#endif  // SAVE_PHOTON_HISTORY
#endif  // STOP_PHOTONS_ON_DETECTION
                                       geoLayerToOMNumIndexPerStringSetLocal);
        }
    }

    // onecellsyncthreads
}

__device__ __forceinline__ void checkForCollision_InCells(const float photonDirLenXYSqr, const float4 &photonPosAndTime,
                                                          const float4 &photonDirAndWlen, const float* getWavelengthBias_dataShared,
#ifdef STOP_PHOTONS_ON_DETECTION
                                                          float &thisStepLength, bool &hitRecorded,
                                                          unsigned short &hitOnString, unsigned short &hitOnDom,
#else
                                                          float thisStepLength, float inv_groupvel,
                                                          float photonTotalPathLength, uint32_t photonNumScatters,
                                                          float distanceTraveledInAbsorptionLengths,
                                                          const float4 &photonStartPosAndTime,
                                                          const float4 &photonStartDirAndWlen,
                                                          const I3CLSimStepCuda &step, uint32_t *hitIndex,
                                                          uint32_t maxHitIndex, I3CLSimPhotonCuda *outputPhotons, 
#ifdef SAVE_PHOTON_HISTORY
                                                          float4 *photonHistory, float4 *currentPhotonHistory,
#endif
#endif
                                                          const unsigned short *geoLayerToOMNumIndexPerStringSetLocal)
{
    // using macros and hard-coded names is
    // not really the best thing to do here..
    // replace with a loop sometime.

#ifdef STOP_PHOTONS_ON_DETECTION
#define DO_CHECK(subdetectorNum)                                                                                 \
    checkForCollision_InCell(photonDirLenXYSqr, photonPosAndTime, photonDirAndWlen, getWavelengthBias_dataShared, thisStepLength, hitRecorded, \
                             hitOnString, hitOnDom, geoLayerToOMNumIndexPerStringSetLocal,                       \
                                                                                                                 \
                            geoCellIndex_##subdetectorNum, GEO_CELL_START_X_##subdetectorNum,  \
                             GEO_CELL_START_Y_##subdetectorNum, GEO_CELL_WIDTH_X_##subdetectorNum,               \
                             GEO_CELL_WIDTH_Y_##subdetectorNum, GEO_CELL_NUM_X_##subdetectorNum,                 \
                             GEO_CELL_NUM_Y_##subdetectorNum);
#else  // STOP_PHOTONS_ON_DETECTION
#ifdef SAVE_PHOTON_HISTORY
#define DO_CHECK(subdetectorNum)                                                                                       \
    checkForCollision_InCell(photonDirLenXYSqr, photonPosAndTime, photonDirAndWlen,getWavelengthBias_dataShared, thisStepLength, inv_groupvel,      \
                             photonTotalPathLength, photonNumScatters, distanceTraveledInAbsorptionLengths,            \
                             photonStartPosAndTime, photonStartDirAndWlen, step, hitIndex, maxHitIndex, outputPhotons \
                             photonHistory, currentPhotonHistory, geoLayerToOMNumIndexPerStringSetLocal,               \
                                                                                                                       \
                             geoCellIndex_##subdetectorNum, GEO_CELL_START_X_##subdetectorNum,                         \
                             GEO_CELL_START_Y_##subdetectorNum, GEO_CELL_WIDTH_X_##subdetectorNum,                     \
                             GEO_CELL_WIDTH_Y_##subdetectorNum, GEO_CELL_NUM_X_##subdetectorNum,                       \
                             GEO_CELL_NUM_Y_##subdetectorNum);
#else  // SAVE_PHOTON_HISTORY
#define DO_CHECK(subdetectorNum)                                                                                    \
    checkForCollision_InCell(                                                                                       \
        photonDirLenXYSqr, photonPosAndTime, photonDirAndWlen,getWavelengthBias_dataShared, thisStepLength, inv_groupvel, photonTotalPathLength, \
        photonNumScatters, distanceTraveledInAbsorptionLengths, photonStartPosAndTime, photonStartDirAndWlen, step, \
        hitIndex, maxHitIndex, outputPhotons, geoLayerToOMNumIndexPerStringSetLocal,                                \
                                                                                                                    \
        geoCellIndex_##subdetectorNum, GEO_CELL_START_X_##subdetectorNum, GEO_CELL_START_Y_##subdetectorNum,        \
        GEO_CELL_WIDTH_X_##subdetectorNum, GEO_CELL_WIDTH_Y_##subdetectorNum, GEO_CELL_NUM_X_##subdetectorNum,      \
        GEO_CELL_NUM_Y_##subdetectorNum);
#endif  // SAVE_PHOTON_HISTORY
#endif  // STOP_PHOTONS_ON_DETECTION

    // argh..
#if GEO_CELL_NUM_SUBDETECTORS > 0
    DO_CHECK(0);
#endif

#if GEO_CELL_NUM_SUBDETECTORS > 1
    DO_CHECK(1);
#endif

#if GEO_CELL_NUM_SUBDETECTORS > 2
    DO_CHECK(2);
#endif

#if GEO_CELL_NUM_SUBDETECTORS > 3
    DO_CHECK(3);
#endif

#if GEO_CELL_NUM_SUBDETECTORS > 4
    DO_CHECK(4);
#endif

#if GEO_CELL_NUM_SUBDETECTORS > 5
    DO_CHECK(5);
#endif

#if GEO_CELL_NUM_SUBDETECTORS > 6
    DO_CHECK(6);
#endif

#if GEO_CELL_NUM_SUBDETECTORS > 7
    DO_CHECK(7);
#endif

#if GEO_CELL_NUM_SUBDETECTORS > 8
    DO_CHECK(8);
#endif

#if GEO_CELL_NUM_SUBDETECTORS > 9
#error more than 9 subdetectors are currently not supported.
#endif

#undef DO_CHECK
}

__device__ __forceinline__ bool checkForCollision(const I3CLPhoton& photon, const I3CLInitialPhoton& photonInitials,
                                                  const I3CLSimStepCuda &step, float &thisStepLength,
                                                  uint32_t *hitIndex, uint32_t maxHitIndex,
                                                  I3CLSimPhotonCuda *outputPhotons,  
                                                  const unsigned short *geoLayerToOMNumIndexPerStringSetLocal, 
                                                  const float* getWavelengthBias_dataShared, int* collisionCheckIndex, CollisionCheck* collisionChecks)
{
    uint32_t myIndex = atomicAdd(collisionCheckIndex, 1);

    // check for collisions
    const float photonDirLenXYSqr = sqr(photon.dirAndWlen.x) + sqr(photon.dirAndWlen.y);
    if (photonDirLenXYSqr <= ZERO)
    {
        collisionChecks[myIndex] = CollisionCheck{ float3{photon.posAndTime.x, photon.posAndTime.y, photon.posAndTime.z},  
                                                    float3{photon.dirAndWlen.x, photon.dirAndWlen.y, photon.dirAndWlen.z},
                                                    thisStepLength, -1, -1};
        return false;
    }


    bool hitRecorded = false;
    unsigned short hitOnString;
    unsigned short hitOnDom;

    float distanceTraveledInAbsorptionLengths = photonInitials.absLength - photon.absLength;
    checkForCollision_InCells(photonDirLenXYSqr, photon.posAndTime, photon.dirAndWlen, getWavelengthBias_dataShared,
                              thisStepLength, hitRecorded, hitOnString, hitOnDom,
                              geoLayerToOMNumIndexPerStringSetLocal);

    // In case photons are stopped on detection
    // (i.e. absorbed by the DOM), we need to record
    // them here (after all possible DOM intersections
    // have been checked).

    if (hitRecorded) {
        saveHit(photon.posAndTime, photon.dirAndWlen, getWavelengthBias_dataShared, thisStepLength, photon.invGroupvel, photon.totalPathLength,
                photon.numScatters, distanceTraveledInAbsorptionLengths, photonInitials.posAndTime, photonInitials.dirAndWlen,
                step, hitOnString, hitOnDom, hitIndex, maxHitIndex, outputPhotons);
    }

    collisionChecks[myIndex] = CollisionCheck{ float3{photon.posAndTime.x, photon.posAndTime.y, photon.posAndTime.z},  
                                                float3{photon.dirAndWlen.x, photon.dirAndWlen.y, photon.dirAndWlen.z},
                                                thisStepLength, (hitRecorded ? hitOnString : -1), (hitRecorded ? hitOnDom : -1)};

    return hitRecorded;
}

#ifdef SAVE_ALL_PHOTONS
#ifdef STOP_PHOTONS_ON_DETECTION
#error The SAVE_ALL_PHOTONS and STOP_PHOTONS_ON_DETECTION options cannot be used at the same time.
#endif
#endif

#ifdef USE_FABS_WORKAROUND
__device__ __forceinline__ float my_fabs(const float a) { return (a < ZERO) ? (-a) : (a); }
#else
__device__ __forceinline__ float my_fabs(const float a) { return fabs(a); }
#endif
__device__ __forceinline__ float sqr(const float a) { return a * a; }

__device__ __forceinline__ int findLayerForGivenZPos(float posZ)
{
    return int((posZ - (float)MEDIUM_LAYER_BOTTOM_POS) / (float)MEDIUM_LAYER_THICKNESS);
}

__device__ __forceinline__ float mediumLayerBoundary(int layer)
{
    return (float(layer) * ((float)MEDIUM_LAYER_THICKNESS)) + (float)MEDIUM_LAYER_BOTTOM_POS;
}

__device__ __forceinline__ void scatterDirectionByAngle(float cosa, float sina, float4 &direction, float randomNumber)
{
    // randomize direction of scattering (rotation around old direction axis)
    const float b = 2.0f * PI * randomNumber;

    const float cosb = cosf(b);
    const float sinb = sinf(b);

    // Rotate new direction into absolute frame of reference
    const float sinth = sqrtf(max(ZERO, ONE - (direction).z * (direction).z));

    if (sinth > 0.f) {  // Current direction not vertical, so rotate
        const float4 oldDir = direction;

        (direction).x = oldDir.x * cosa - ((oldDir.y * cosb + oldDir.z * oldDir.x * sinb) * sina) / sinth;
        (direction).y = oldDir.y * cosa + ((oldDir.x * cosb - oldDir.z * oldDir.y * sinb) * sina) / sinth;
        (direction).z = oldDir.z * cosa + sina * sinb * sinth;
    } else {  // Current direction is vertical, so this is trivial
        (direction).x = sina * cosb;
        (direction).y = sina * sinb;
        (direction).z = cosa * sign((direction).z);
    }

    {
        const float recip_length = rsqrtf(sqr((direction).x) + sqr((direction).y) + sqr((direction).z));

        (direction).x *= recip_length;
        (direction).y *= recip_length;
        (direction).z *= recip_length;
    }

    // printf("direction after=(%f,%f,%f) len^2=%f\n",
    //       (*direction).x, (*direction).y, (*direction).z,
    //       (*direction).x*(*direction).x + (*direction).y*(*direction).y +
    //       (*direction).z*(*direction).z);
}

__device__ __forceinline__ void createPhotonFromTrack(const I3CLSimStepCuda &step, const float4 &stepDir, RNG_ARGS,
                                                      float4 &photonPosAndTime, float4 &photonDirAndWlen,
                                                      float* _generateWavelength_0distYValuesShared, float* _generateWavelength_0distYCumulativeValuesShared)
{
    float shiftMultiplied = step.dirAndLengthAndBeta.z * RNG_CALL_UNIFORM_CO;
    float inverseParticleSpeed = 1.f / (speedOfLight * step.dirAndLengthAndBeta.w);

    // move along the step direction
    photonPosAndTime = float4{
        step.posAndTime.x + stepDir.x * shiftMultiplied, step.posAndTime.y + stepDir.y * shiftMultiplied,
        step.posAndTime.z + stepDir.z * shiftMultiplied, step.posAndTime.w + inverseParticleSpeed * shiftMultiplied};

    // determine the photon layer (clamp if necessary)
    unsigned int layer = min(max(findLayerForGivenZPos((photonPosAndTime).z), 0), MEDIUM_LAYERS - 1);

#ifndef NO_FLASHER
    if (step.sourceType == 0) {
#endif
        // sourceType==0 is always Cherenkov light with the correct angle w.r.t. the
        // particle/step

        // our photon still needs a wavelength. create one!
        const float wavelength = generateWavelength_0(RNG_ARGS_TO_CALL, _generateWavelength_0distYValuesShared,   _generateWavelength_0distYCumulativeValuesShared);

        const float cosCherenkov = min(
            ONE, 1.f / (step.dirAndLengthAndBeta.w * getPhaseRefIndex(layer, wavelength)));  // cos theta = 1/(beta*n)
        const float sinCherenkov = sqrtf(ONE - cosCherenkov * cosCherenkov);
        // determine the photon direction

        // start with the track direction
        (photonDirAndWlen).x = stepDir.x;
        (photonDirAndWlen).y = stepDir.y;
        (photonDirAndWlen).z = stepDir.z;
        (photonDirAndWlen).w = wavelength;

        // and now rotate to cherenkov emission direction
        scatterDirectionByAngle(cosCherenkov, sinCherenkov, photonDirAndWlen, RNG_CALL_UNIFORM_CO);

#ifndef NO_FLASHER
    } else {
        // steps >= 1 are flasher emissions, they do not need cherenkov rotation

        const float wavelength = generateWavelength(uint(step.sourceType), RNG_ARGS_TO_CALL, _generateWavelength_0distYValuesShared,   _generateWavelength_0distYCumulativeValuesShared);

        // use the step direction as the photon direction
        (photonDirAndWlen).x = stepDir.x;
        (photonDirAndWlen).y = stepDir.y;
        (photonDirAndWlen).z = stepDir.z;
        (photonDirAndWlen).w = wavelength;
    }
#endif
}

__device__ __forceinline__ float2 sphDirFromCar(float4 carDir)
{
    // Calculate Spherical coordinates from Cartesian
    const float r_inv = rsqrtf(carDir.x * carDir.x + carDir.y * carDir.y + carDir.z * carDir.z);

    float theta = 0.f;
    if ((my_fabs(carDir.z * r_inv)) <= 1.f) {
        theta = acos(carDir.z * r_inv);
    } else {
        if (carDir.z < 0.f) theta = CUDART_PI_F;
    }
    if (theta < 0.f) theta += 2.f * CUDART_PI_F;

    float phi = atan2(carDir.y, carDir.x);
    if (phi < 0.f) phi += 2.f * CUDART_PI_F;

    return float2{theta, phi};
}

// Record a photon on a DOM
__device__ __forceinline__ void saveHit(const float4 &photonPosAndTime, const float4 &photonDirAndWlen, const float* getWavelengthBias_dataShared,
                                        const float thisStepLength, float inv_groupvel, float photonTotalPathLength,
                                        uint32_t photonNumScatters, float distanceTraveledInAbsorptionLengths,
                                        const float4 &photonStartPosAndTime, const float4 &photonStartDirAndWlen,
                                        const I3CLSimStepCuda &step, unsigned short hitOnString,
                                        unsigned short hitOnDom,
                                        uint32_t *hitIndex,  // shared
                                        uint32_t maxHitIndex, I3CLSimPhotonCuda *outputPhotons                                 
#ifdef SAVE_PHOTON_HISTORY
                                        ,
                                        float4 *photonHistory, float4 *currentPhotonHistory
#endif
)
{
    // PRINTL
    uint32_t myIndex = atomicAdd(&hitIndex[0], 1);

    if (myIndex < maxHitIndex) {
        // Emit photon position relative to the hit DOM
#ifndef CABLE_RADIUS
        float domPosX, domPosY, domPosZ;
        geometryGetDomPosition(hitOnString, hitOnDom, domPosX, domPosY, domPosZ);
#else
        float domPosX, domPosY, domPosZ, cableOrientation;
        geometryGetDomPosition(hitOnString, hitOnDom, domPosX, domPosY, domPosZ, &cableOrientation);
#endif
#ifdef PANCAKE_FACTOR
        {
            // undo pancaking by scaling the distance of closest approach to the
            // DOM center
            float px = photonPosAndTime.x - domPosX;
            float py = photonPosAndTime.y - domPosY;
            float pz = photonPosAndTime.z - domPosZ;
            float parallel = px * photonDirAndWlen.x + py * photonDirAndWlen.y + pz * photonDirAndWlen.z;
            float nx = px - parallel * photonDirAndWlen.x;
            float ny = py - parallel * photonDirAndWlen.y;
            float nz = pz - parallel * photonDirAndWlen.z;
            domPosX += ((PANCAKE_FACTOR - ONE) / PANCAKE_FACTOR) * nx;
            domPosY += ((PANCAKE_FACTOR - ONE) / PANCAKE_FACTOR) * ny;
            domPosZ += ((PANCAKE_FACTOR - ONE) / PANCAKE_FACTOR) * nz;
        }
#endif
        I3CLSimPhotonCuda outphoton;
        outphoton.posAndTime = float4{photonPosAndTime.x + thisStepLength * photonDirAndWlen.x - domPosX,
                                      photonPosAndTime.y + thisStepLength * photonDirAndWlen.y - domPosY,
                                      photonPosAndTime.z + thisStepLength * photonDirAndWlen.z - domPosZ,
                                      photonPosAndTime.w + thisStepLength * inv_groupvel};

        outphoton.dir = sphDirFromCar(photonDirAndWlen);
        outphoton.wavelength = photonDirAndWlen.w;

        outphoton.cherenkovDist = photonTotalPathLength + thisStepLength;
        outphoton.numScatters = photonNumScatters;
        outphoton.weight = step.weight / getWavelengthBias(photonDirAndWlen.w, getWavelengthBias_dataShared);
        outphoton.identifier = step.identifier;

        outphoton.stringID = short(hitOnString);
        outphoton.omID = ushort(hitOnDom);

        outphoton.startPosAndTime = photonStartPosAndTime;

        outphoton.startDir = sphDirFromCar(photonStartDirAndWlen);

        outphoton.groupVelocity = 1.f / (inv_groupvel);

        outphoton.distInAbsLens = distanceTraveledInAbsorptionLengths;

        outputPhotons[myIndex] = outphoton;

#ifdef SAVE_PHOTON_HISTORY
        for (uint32_t i = 0; i < NUM_PHOTONS_IN_HISTORY; ++i) {
            photonHistory[myIndex * NUM_PHOTONS_IN_HISTORY + i] = currentPhotonHistory[i];
        }
#endif
    }
}

#endif  // PROPAGATIONKERNELFUNCTIUONS_CUH