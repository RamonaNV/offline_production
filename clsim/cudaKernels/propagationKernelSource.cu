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

/* 
    implements main simulation kernel as well as host code to launch it
*/

// includes
// ------------------
#include "propagationKernelSource.cuh"

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cooperative_groups.h>
#include <cooperative_groups/memcpy_async.h>
namespace cg = cooperative_groups;

#include <iostream>
#include <chrono>

#include "settings.cuh"
#include "dataStructCuda.cuh"
#include "utils.cuh"
#include "rng.cuh"
#include "propagationKernelFunctions.cuh"
#include "zOffsetHandling.cuh"
#include "wlenGeneration.cuh"
#include "scatteringAndAbsorbtionData.cuh"
// ------------------

// remark: ignored tabulate version, removed ifdef TABULATE
// also removed ifdef DOUBLEPRECISION.
// SAVE_PHOTON_HISTORY  and SAVE_ALL_PHOTONS are not define for now, i.e. commented out these snippets,
// s.t. it corresponds to the default contstructor of I3CLSimStepToPhotonConverterOpenCL

__global__ __launch_bounds__(NTHREADS_PER_BLOCK, NBLOCKS_PER_SM) void propKernel( I3CLSimStepCuda* __restrict__ steps, int numSteps, 
                                                                    uint32_t* hitIndex, uint32_t maxHitIndex, I3CLSimPhotonCuda* __restrict__ outputPhotons,
                                                                    const float* wlenLut, const float* zOffsetLut, 
                                                                    const unsigned short* __restrict__ geoLayerToOMNumIndexPerStringSet, 
                                                                    uint64_t* __restrict__ rng_x, uint32_t* __restrict__ rng_a); 

void launch_CudaPropogate(const I3CLSimStep* __restrict__ in_steps, int nsteps, const uint32_t maxHitIndex,
                          unsigned short* geoLayerToOMNumIndexPerStringSet, int ngeolayer,
                          I3CLSimPhotonSeries& outphotons, uint64_t* __restrict__ MWC_RNG_x,
                          uint32_t* __restrict__ MWC_RNG_a, int sizeRNG, float& totalCudaKernelTime)
{
    // setup the rng
    uint64_t* d_MWC_RNG_x;
    uint32_t* d_MWC_RNG_a;
    initMWCRng(sizeRNG, MWC_RNG_x, MWC_RNG_a, &d_MWC_RNG_x, &d_MWC_RNG_a);

    printf("nsteps total = %d but dividing into %d launches of max size %d \n", nsteps, 1, nsteps);

    // upload "geo layer per string set" data 
    unsigned short* d_geolayer;
    CUDA_ERR_CHECK(cudaMalloc((void**)&d_geolayer, ngeolayer * sizeof(unsigned short)));
    CUDA_ERR_CHECK(cudaMemcpy(d_geolayer, geoLayerToOMNumIndexPerStringSet, ngeolayer * sizeof(unsigned short),
                              cudaMemcpyHostToDevice));

    // convert and upload steps
    I3CLSimStepCuda* h_cudastep = (I3CLSimStepCuda*)malloc(nsteps * sizeof(struct I3CLSimStepCuda));
    for (int i = 0; i < nsteps; i++) {
        h_cudastep[i] = I3CLSimStepCuda(in_steps[i]);
    }
    I3CLSimStepCuda* d_cudastep;
    CUDA_ERR_CHECK(cudaMalloc((void**)&d_cudastep, nsteps * sizeof(I3CLSimStepCuda)));
    CUDA_ERR_CHECK(cudaMemcpy(d_cudastep, h_cudastep, nsteps * sizeof(I3CLSimStepCuda), cudaMemcpyHostToDevice));

    // allocate storage to store hits
    uint32_t* d_hitIndex;
    uint32_t h_hitIndex[1];
    h_hitIndex[0] = 0;
    CUDA_ERR_CHECK(cudaMalloc((void**)&d_hitIndex, 1 * sizeof(uint32_t)));
    CUDA_ERR_CHECK(cudaMemcpy(d_hitIndex, h_hitIndex, 1 * sizeof(uint32_t), cudaMemcpyHostToDevice));

    I3CLSimPhotonCuda* d_cudaphotons;
    CUDA_ERR_CHECK(cudaMalloc((void**)&d_cudaphotons, maxHitIndex * sizeof(I3CLSimPhotonCuda)));

    // wlen lut
    auto wlenLut = generateWavelengthLut();
    float* d_wlenLut;
    CUDA_ERR_CHECK(cudaMalloc((void**)&d_wlenLut, WLEN_LUT_SIZE * sizeof(float)));
    CUDA_ERR_CHECK(cudaMemcpy(d_wlenLut, wlenLut.data(), WLEN_LUT_SIZE * sizeof(float), cudaMemcpyHostToDevice));

    // zOffset lut
    auto zOffsetLut = generateZOffsetLut();
    float* d_zOffsetLut;
    CUDA_ERR_CHECK(cudaMalloc((void**)&d_zOffsetLut, zOffsetLut.size() * sizeof(float)));
    CUDA_ERR_CHECK(cudaMemcpy(d_zOffsetLut, zOffsetLut.data(), zOffsetLut.size() * sizeof(float), cudaMemcpyHostToDevice));

    // compute block number and launch
    int numBlocks = (nsteps + NTHREADS_PER_BLOCK - 1) / NTHREADS_PER_BLOCK;
    printf("launching kernel propKernel<<< %d , %d >>>( .., nsteps=%d)  \n", numBlocks, NTHREADS_PER_BLOCK, nsteps);

    std::chrono::time_point<std::chrono::system_clock> startKernel = std::chrono::system_clock::now();
    propKernel<<<numBlocks, NTHREADS_PER_BLOCK>>>(d_cudastep, nsteps, 
                                                  d_hitIndex, maxHitIndex, d_cudaphotons,
                                                  d_wlenLut, d_zOffsetLut,
                                                  d_geolayer,
                                                  d_MWC_RNG_x, d_MWC_RNG_a);

    CUDA_ERR_CHECK(cudaDeviceSynchronize());
    std::chrono::time_point<std::chrono::system_clock> endKernel = std::chrono::system_clock::now();
    totalCudaKernelTime = std::chrono::duration_cast<std::chrono::milliseconds>(endKernel - startKernel).count();

    CUDA_ERR_CHECK(cudaMemcpy(h_hitIndex, d_hitIndex, 1 * sizeof(uint32_t), cudaMemcpyDeviceToHost));
    int numberPhotons = h_hitIndex[0];

    if (numberPhotons > maxHitIndex) {
        printf("Maximum number of photons exceeded, only receiving %u of %u photons", maxHitIndex, numberPhotons);
        numberPhotons = maxHitIndex;
    }

    // copy (max fo maxHitIndex) photons to host.
    struct I3CLSimPhotonCuda* h_cudaphotons =
        (struct I3CLSimPhotonCuda*)malloc(numberPhotons * sizeof(struct I3CLSimPhotonCuda));
    CUDA_ERR_CHECK(
        cudaMemcpy(h_cudaphotons, d_cudaphotons, numberPhotons * sizeof(I3CLSimPhotonCuda), cudaMemcpyDeviceToHost));

    outphotons.resize(numberPhotons);
    for (int i = 0; i < numberPhotons; i++) {
        outphotons[i] = h_cudaphotons[i].getI3CLSimPhoton();
    }

    free(h_cudastep);
    free(h_cudaphotons);
    cudaFree(d_cudaphotons);
    cudaFree(d_cudastep);
    cudaFree(d_geolayer);
    cudaFree(d_MWC_RNG_a);
    cudaFree(d_MWC_RNG_x);
    printf("photon hits = %i from %i steps \n", numberPhotons, nsteps);
}

__global__ void propKernel( I3CLSimStepCuda* __restrict__ steps, int numSteps, 
                            uint32_t* hitIndex, uint32_t maxHitIndex, I3CLSimPhotonCuda* __restrict__ outputPhotons,
                            const float* wlenLut, const float* zOffsetLut, 
                            const unsigned short* __restrict__ geoLayerToOMNumIndexPerStringSet, 
                            uint64_t* __restrict__ rng_x, uint32_t* __restrict__ rng_a)
{
     // copy some LUTs to shared memory for faster access
    cg::thread_block block = cg::this_thread_block();

    #ifdef SHARED_WLEN
        __shared__ float sharedWlenLut[WLEN_LUT_SIZE];
        cg::memcpy_async(block, sharedWlenLut, wlenLut, sizeof(float)*WLEN_LUT_SIZE);
        const float* wlenLutPointer = sharedWlenLut;
    #else
        const float* wlenLutPointer = wlenLut;
    #endif

    #ifdef SHARED_ZOFFSET
        __shared__ float sharedZOffsetLut[ZOLUT_SIZE];
        cg::memcpy_async(block, sharedZOffsetLut, zOffsetLut, sizeof(float)*ZOLUT_SIZE);
        const float* zOffsetLutPointer = sharedZOffsetLut;
    #else
        const float* zOffsetLutPointer = zOffsetLut;
    #endif

    #ifdef SHARED_ICE_PROPERTIES
        __shared__ float sharedScatteringLength[171];
        __shared__ float sharedAbsorptionADust[171];
        __shared__ float sharedAbsorptionDeltaTau[171];
        cg::memcpy_async(block, sharedScatteringLength, scatteringLength_b400_LUT, sizeof(float)*171);
        cg::memcpy_async(block, sharedAbsorptionADust, absorptionLength_aDust400_LUT, sizeof(float)*171);
        cg::memcpy_async(block, sharedAbsorptionDeltaTau, absorptionLength_deltaTau_LUT, sizeof(float)*171);
        const float* scatteringLutPointer = sharedScatteringLength;
        const float* absorbtionLutPointer = sharedAbsorptionADust;
        const float* absorbtionDeltaTauLutPointer = sharedAbsorptionDeltaTau;
    #else
        const float* scatteringLutPointer = scatteringLength_b400_LUT;
        const float* absorbtionLutPointer = absorptionLength_aDust400_LUT;
        const float* absorbtionDeltaTauLutPointer = absorptionLength_deltaTau_LUT;
    #endif

    #ifdef SHARED_WLEN_BIAS
        __shared__ float getWavelengthBias_dataShared[43];
        cg::memcpy_async(block, getWavelengthBias_dataShared, getWavelengthBias_data, sizeof(float)*43);
        const float* wlenBiasLutPointer = getWavelengthBias_dataShared;
    #else
        const float* wlenBiasLutPointer = getWavelengthBias_data;
    #endif

    #ifdef SHARED_NUM_INDEX_STRING_SET
        __shared__ unsigned short geoLayerToOMNumIndexPerStringSetLocal[GEO_geoLayerToOMNumIndexPerStringSet_BUFFER_SIZE];
        cg::memcpy_async(block, geoLayerToOMNumIndexPerStringSetLocal, geoLayerToOMNumIndexPerStringSet, sizeof(unsigned short)*GEO_geoLayerToOMNumIndexPerStringSet_BUFFER_SIZE);
        const unsigned short* numIndexStringSetPointer = geoLayerToOMNumIndexPerStringSetLocal;
    #else
        const unsigned short* numIndexStringSetPointer = geoLayerToOMNumIndexPerStringSet;
    #endif

    #ifdef SHARED_COLLISION_GRID_DATA
        __shared__ unsigned short geoCellIndex0shared[GEO_CELL_NUM_X_0 * GEO_CELL_NUM_Y_0];
        __shared__ unsigned short geoCellIndex1shared[GEO_CELL_NUM_X_1 * GEO_CELL_NUM_Y_1];
        cg::memcpy_async(block, geoCellIndex0shared, geoCellIndex_0, sizeof(unsigned short)*GEO_CELL_NUM_X_0 * GEO_CELL_NUM_Y_0);
        cg::memcpy_async(block, geoCellIndex1shared, geoCellIndex_1, sizeof(unsigned short)*GEO_CELL_NUM_X_1 * GEO_CELL_NUM_Y_1);
        const unsigned short* geoCellIndex0Pointer = geoCellIndex0shared;
        const unsigned short* geoCellIndex1Pointer = geoCellIndex1shared;
    #else
        const unsigned short* geoCellIndex0Pointer = geoCellIndex_0;
        const unsigned short* geoCellIndex1Pointer = geoCellIndex_1;
    #endif

    #ifdef SHARED_STRING_DATA
        __shared__ unsigned char geoStringInSetShared[NUM_STRINGS];
        __shared__ unsigned short geoLayerNumShared[GEO_LAYER_STRINGSET_NUM];
        __shared__ float geoLayerStartZShared[GEO_LAYER_STRINGSET_NUM];
        __shared__ float geoLayerHeightShared[GEO_LAYER_STRINGSET_NUM];
        cg::memcpy_async(block, geoStringInSetShared, geoStringInStringSet, sizeof(unsigned char)*NUM_STRINGS);
        cg::memcpy_async(block, geoLayerNumShared, geoLayerNum, sizeof(unsigned short)*GEO_LAYER_STRINGSET_NUM);
        cg::memcpy_async(block, geoLayerStartZShared, geoLayerStartZ, sizeof(float)*GEO_LAYER_STRINGSET_NUM);
        cg::memcpy_async(block, geoLayerHeightShared, geoLayerHeight, sizeof(float)*GEO_LAYER_STRINGSET_NUM);
        const unsigned char* geoStringInSetPointer = geoStringInSetShared;
        const unsigned short* geoLayerNumPointer = geoLayerNumShared;
        const float* geoLayerStartZPointer = geoLayerStartZShared;
        const float* geoLayerHeightPointer = geoLayerHeightShared;
    #else
        const unsigned char* geoStringInSetPointer = geoStringInStringSet;
        const unsigned short* geoLayerNumPointer = geoLayerNum;
        const float* geoLayerStartZPointer = geoLayerStartZ;
        const float* geoLayerHeightPointer = geoLayerHeight;
    #endif

    cg::wait(block);

    // get thread id
    int id = threadIdx.x + blockIdx.x * blockDim.x;
    if(id > numSteps)
        return;

    // initialize rng
    RngType rng(rng_x[id],rng_a[id]);

    // load step and calculate direction
    const I3CLSimStepCuda step = steps[id];
    const float3 stepDir = calculateStepDir(step);

    // variables to store data about current photon
    uint32_t photonsLeftToPropagate = step.numPhotons;
    I3CLPhoton photon;
    photon.absLength = 0.0f;

    // loop until all photons are done
    while (photonsLeftToPropagate > 0) {

        // if current photon is done, create a new one
        if (photon.absLength < EPSILON) {
            photon = createPhoton(step, stepDir, wlenLutPointer, rng);
        }

        // propagate through layers until scattered or absorbed
        float distanceTraveled;
        bool absorbed = propPhoton(photon, distanceTraveled, rng, scatteringLutPointer, absorbtionLutPointer, absorbtionDeltaTauLutPointer, zOffsetLutPointer);

        // check for collision with DOMs, if collision has happened, the hit will be stored in outputPhotons
        bool collided = checkForCollisionOld(photon, step, distanceTraveled, 
                                  hitIndex, maxHitIndex, outputPhotons, numIndexStringSetPointer, wlenBiasLutPointer, 
                                  geoCellIndex0Pointer, geoCellIndex1Pointer, geoStringInSetPointer, geoLayerNumPointer, geoLayerStartZPointer, geoLayerHeightPointer);

        // remove photon if it is collided or absorbed
        // we get the next photon at the beginning of the loop
        if (collided || absorbed) {
            photon.absLength = 0.0f;
            --photonsLeftToPropagate;
        }
        else
        {
            // move the photon along its current direction for the distance it was propagated through the ice
            // then scatter to find a new direction vector
            updatePhotonTrack(photon, distanceTraveled);
            scatterPhoton(photon, rng);
        }
    }

    // store rng state
    rng_x[id] = rng.getState();
}
