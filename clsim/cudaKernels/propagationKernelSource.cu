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
    implements main simulation kernel as well as host code to launch it
*/

// includes
// ------------------
#include "propagationKernelSource.cuh"

#include <cuda.h>
#include <cuda_runtime_api.h>
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
                                                                    uint64_t* __restrict__ rng_x, uint32_t* __restrict__ rng_a); 

__global__ __launch_bounds__(NTHREADS_PER_BLOCK, NBLOCKS_PER_SM) void propKernelJobqueue(I3CLSimStepCuda* __restrict__ steps, int numSteps, 
                                                                    uint32_t* hitIndex, uint32_t maxHitIndex, I3CLSimPhotonCuda* __restrict__ outputPhotons,
                                                                    const float* wlenLut, const float* zOffsetLut,
                                                                    uint64_t* __restrict__ rng_x, uint32_t* __restrict__ rng_a, int numPrimes); 

struct KernelBuffers {
    float* wlenLut;
    float* zOffsetLut;
    uint64_t *MWC_RNG_x;
    uint32_t *MWC_RNG_a;
    I3CLSimStepCuda* inputSteps;
    uint32_t numInputSteps;
    I3CLSimPhotonCuda* outputPhotons;
    uint32_t *numOutputPhotons;
    uint32_t maxHitIndex;
    cudaStream_t stream;
    KernelBuffers(
        size_t maxNumWorkItems, size_t maxNumOutputPhotons,
        const std::vector<uint64_t> &x, const std::vector<uint32_t> &a
    ) : wlenLut(nullptr), zOffsetLut(nullptr), MWC_RNG_x(nullptr), MWC_RNG_a(nullptr), inputSteps(nullptr), numInputSteps(0), outputPhotons(nullptr), numOutputPhotons(nullptr), maxHitIndex(maxNumOutputPhotons) {
        {
            auto wlen = generateWavelengthLut();
            CUDA_ERR_THROW(cudaMalloc((void**)&wlenLut, wlen.size()*sizeof(float)));
            CUDA_ERR_THROW(cudaMemcpy(wlenLut, wlen.data(), wlen.size() * sizeof(float), cudaMemcpyHostToDevice));
        }
        {
            auto zOffset = generateZOffsetLut();
            CUDA_ERR_THROW(cudaMalloc((void**)&zOffsetLut, zOffset.size() * sizeof(float)));
            CUDA_ERR_THROW(cudaMemcpy(zOffsetLut, zOffset.data(), zOffset.size() * sizeof(float), cudaMemcpyHostToDevice));
        }
        {
            // FIXME: add jobqueue support back
            initMWCRng(x.size(), x.data(), a.data(), &MWC_RNG_x, &MWC_RNG_a);
        }
        CUDA_ERR_THROW(cudaMalloc((void**)&inputSteps, maxNumWorkItems * sizeof(I3CLSimStepCuda)));
        CUDA_ERR_THROW(cudaMalloc((void**)&outputPhotons, maxNumOutputPhotons * sizeof(I3CLSimPhotonCuda)));
        CUDA_ERR_THROW(cudaMalloc((void**)&numOutputPhotons, sizeof(uint32_t)));
        CUDA_ERR_THROW(cudaStreamCreate(&stream));
    }

    ~KernelBuffers() {
        if (wlenLut != nullptr)
            cudaFree(wlenLut);
        if (zOffsetLut != nullptr)
            cudaFree(zOffsetLut);
        if (MWC_RNG_x != nullptr)
            cudaFree(MWC_RNG_x);
        if (MWC_RNG_a != nullptr)
            cudaFree(MWC_RNG_a);
        if (inputSteps != nullptr)
            cudaFree(inputSteps);
        if (outputPhotons != nullptr)
            cudaFree(outputPhotons);
        if (numOutputPhotons != nullptr)
            cudaFree(numOutputPhotons);
        cudaStreamDestroy(stream);
    }

// hide from device compilation trajectory (I3CLSimStep contains unsupported vector types)
#ifndef __CUDA_ARCH__
    void uploadSteps(const std::vector<I3CLSimStep> &steps) {
        std::vector<I3CLSimStepCuda> cudaSteps(steps.size());
        for (int i = 0; i < steps.size(); i++) {
            cudaSteps[i] = I3CLSimStepCuda(steps[i]);
        }
        CUDA_ERR_THROW(cudaMemcpyAsync(inputSteps, cudaSteps.data(), cudaSteps.size() * sizeof(I3CLSimStepCuda), cudaMemcpyHostToDevice, stream));
        CUDA_ERR_THROW(cudaStreamSynchronize(stream));
        numInputSteps = cudaSteps.size();
    }
    std::vector<I3CLSimPhoton> downloadPhotons() {
        uint32_t numberPhotons;
        CUDA_ERR_THROW(cudaMemcpyAsync(&numberPhotons, numOutputPhotons, 1 * sizeof(uint32_t), cudaMemcpyDeviceToHost, stream));
        CUDA_ERR_THROW(cudaStreamSynchronize(stream));
        std::vector<I3CLSimPhotonCuda> cudaPhotons(numberPhotons);
        std::vector<I3CLSimPhoton> photons(numberPhotons);
        CUDA_ERR_THROW(cudaMemcpyAsync(cudaPhotons.data(), outputPhotons, numberPhotons * sizeof(I3CLSimPhotonCuda), cudaMemcpyDeviceToHost, stream));
        CUDA_ERR_THROW(cudaStreamSynchronize(stream));
        for (int i = 0; i < numberPhotons; i++) {
            photons[i] = cudaPhotons[i].getI3CLSimPhoton();
        }
        return photons;
    }
#endif
};

Kernel::Kernel(
    int device,
    size_t maxNumWorkItems,
    size_t maxNumOutputPhotons,
    const std::vector<uint64_t> &x,
    const std::vector<uint32_t> &a
) : impl(new KernelBuffers(maxNumWorkItems, maxNumOutputPhotons, x, a))
{
    CUDA_ERR_THROW(cudaSetDevice(device));
}

// dtor here, where KernelBuffers is complete
Kernel::~Kernel() {}

#ifndef __CUDA_ARCH__
void Kernel::uploadSteps(const std::vector<I3CLSimStep> &steps) { impl->uploadSteps(steps); }
std::vector<I3CLSimPhoton> Kernel::downloadPhotons() { return impl->downloadPhotons(); }
void Kernel::execute() {
#ifdef USE_JOBQUEUE
    propKernelJobqueue<<<numBlocks, NTHREADS_PER_BLOCK, 0, impl->stream >>>(impl->inputSteps, impl->numInputSteps,
                                                impl->numOutputPhotons, impl->maxHitIndex, impl->outputPhotons,
                                                impl->wlenLut, impl->zOffsetLut,
                                                impl->MWC_RNG_x, impl->MWC_RNG_a, sizeRNG);
#else
    int numBlocks = (impl->numInputSteps + NTHREADS_PER_BLOCK - 1) / NTHREADS_PER_BLOCK;
    propKernel<<<numBlocks, NTHREADS_PER_BLOCK, 0, impl->stream >>>(impl->inputSteps, impl->numInputSteps,
                                                impl->numOutputPhotons, impl->maxHitIndex, impl->outputPhotons,
                                                impl->wlenLut, impl->zOffsetLut,
                                                impl->MWC_RNG_x, impl->MWC_RNG_a);
#endif
    CUDA_ERR_THROW(cudaStreamSynchronize(impl->stream));
}
#endif

void launch_CudaPropogate(const I3CLSimStep* __restrict__ in_steps, int nsteps, const uint32_t maxHitIndex,
                          I3CLSimPhotonSeries& outphotons, uint64_t* __restrict__ MWC_RNG_x,
                          uint32_t* __restrict__ MWC_RNG_a, int sizeRNG, float& totalCudaKernelTime)
{
    // setup the rng
    uint64_t* d_MWC_RNG_x;
    uint32_t* d_MWC_RNG_a;
    #ifdef USE_JOBQUEUE
        int numBlocks = (nsteps*32 + NTHREADS_PER_BLOCK - 1) / NTHREADS_PER_BLOCK; // run 32 threads per step
        int numThreads = numBlocks * NTHREADS_PER_BLOCK;
        initMWCRng_jobqueue(sizeRNG, numThreads, MWC_RNG_x, MWC_RNG_a, &d_MWC_RNG_x, &d_MWC_RNG_a);
    #else
         int numBlocks = (nsteps + NTHREADS_PER_BLOCK - 1) / NTHREADS_PER_BLOCK;
        initMWCRng(sizeRNG, MWC_RNG_x, MWC_RNG_a, &d_MWC_RNG_x, &d_MWC_RNG_a);
    #endif

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

    // launch
    printf("launching kernel propKernel<<< %d , %d >>>( .., nsteps=%d)  \n", numBlocks, NTHREADS_PER_BLOCK, nsteps);
    std::chrono::time_point<std::chrono::system_clock> startKernel = std::chrono::system_clock::now();

    #ifdef USE_JOBQUEUE
        propKernelJobqueue<<<numBlocks, NTHREADS_PER_BLOCK>>>(d_cudastep, nsteps,
                                                    d_hitIndex, maxHitIndex, d_cudaphotons,
                                                    d_wlenLut, d_zOffsetLut,
                                                    d_MWC_RNG_x, d_MWC_RNG_a, sizeRNG);
    #else
        propKernel<<<numBlocks, NTHREADS_PER_BLOCK>>>(d_cudastep, nsteps,
                                                    d_hitIndex, maxHitIndex, d_cudaphotons,
                                                    d_wlenLut, d_zOffsetLut,
                                                    d_MWC_RNG_x, d_MWC_RNG_a);
    #endif

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
    CUDA_ERR_CHECK(cudaFree(d_cudaphotons));
    CUDA_ERR_CHECK(cudaFree(d_cudastep));
    CUDA_ERR_CHECK(cudaFree(d_MWC_RNG_a));
    CUDA_ERR_CHECK(cudaFree(d_MWC_RNG_x));
    printf("photon hits = %i from %i steps \n", numberPhotons, nsteps);
}

__global__ void propKernel( I3CLSimStepCuda* __restrict__ steps, int numSteps, 
                            uint32_t* hitIndex, uint32_t maxHitIndex, I3CLSimPhotonCuda* __restrict__ outputPhotons,
                            const float* wlenLut, const float* zOffsetLut,
                            uint64_t* __restrict__ rng_x, uint32_t* __restrict__ rng_a)
{
    #ifdef SHARED_WLEN
        __shared__ float sharedWlenLut[WLEN_LUT_SIZE];
        for (int i = threadIdx.x; i < WLEN_LUT_SIZE; i += blockDim.x) {
            sharedWlenLut[i] = wlenLut[i];
        }
        const float* wlenLutPointer = sharedWlenLut;
    #else
        const float* wlenLutPointer = wlenLut;
    #endif

    #ifdef SHARED_ICE_PROPERTIES
        __shared__ float sharedScatteringLength[171];
        __shared__ float sharedAbsorptionADust[171];
        __shared__ float sharedAbsorptionDeltaTau[171];
        for (int i = threadIdx.x; i < 171; i += blockDim.x) {
            sharedScatteringLength[i] = scatteringLength_b400_LUT[i];
            sharedAbsorptionADust[i] = absorptionLength_aDust400_LUT[i];
            sharedAbsorptionDeltaTau[i] = absorptionLength_deltaTau_LUT[i];
        }
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
        for (int i = threadIdx.x; i < 43; i += blockDim.x) {
            getWavelengthBias_dataShared[i] = getWavelengthBias_data[i];
        }
        const float* wlenBiasLutPointer = getWavelengthBias_dataShared;
    #else
        const float* wlenBiasLutPointer = getWavelengthBias_data;
    #endif

    #ifdef SHARED_NUM_INDEX_STRING_SET
        __shared__ unsigned short geoLayerToOMNumIndexPerStringSetLocal[GEO_geoLayerToOMNumIndexPerStringSet_BUFFER_SIZE];
        for (int i = threadIdx.x; i < GEO_geoLayerToOMNumIndexPerStringSet_BUFFER_SIZE; i += blockDim.x) {
            geoLayerToOMNumIndexPerStringSetLocal[i] = geoLayerToOMNumIndexPerStringSet[i];
        }
        const unsigned short* numIndexStringSetPointer = geoLayerToOMNumIndexPerStringSetLocal;
    #else
        const unsigned short* numIndexStringSetPointer = geoLayerToOMNumIndexPerStringSet;
    #endif

    #ifdef SHARED_COLLISION_GRID_DATA
        __shared__ unsigned short geoCellIndex0shared[GEO_CELL_NUM_X_0 * GEO_CELL_NUM_Y_0];
        __shared__ unsigned short geoCellIndex1shared[GEO_CELL_NUM_X_1 * GEO_CELL_NUM_Y_1];
        for (int i = threadIdx.x; i < GEO_CELL_NUM_X_0 * GEO_CELL_NUM_Y_0; i += blockDim.x) {
            geoCellIndex0shared[i] = geoCellIndex_0[i];
        }
        for (int i = threadIdx.x; i < GEO_CELL_NUM_X_1 * GEO_CELL_NUM_Y_1; i += blockDim.x) {
            geoCellIndex1shared[i] = geoCellIndex_1[i];
        }
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
        for (int i = threadIdx.x; i < NUM_STRINGS; i += blockDim.x) {
            geoStringInSetShared[i] = geoStringInStringSet[i];
        }
        for (int i = threadIdx.x; i < GEO_LAYER_STRINGSET_NUM; i += blockDim.x) {
            geoLayerNumShared[i] = geoLayerNum[i];
            geoLayerStartZShared[i] = geoLayerStartZ[i];
            geoLayerHeightShared[i] = geoLayerHeight[i];
        }
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

    #ifdef SHARED_STRING_POSITIONS
        __shared__ float geoStringPosXShared[NUM_STRINGS];
        __shared__ float geoStringPosYShared[NUM_STRINGS];
        __shared__ float geoStringMinZShared[NUM_STRINGS];
        __shared__ float geoStringMaxZShared[NUM_STRINGS];
        for (int i = threadIdx.x; i < NUM_STRINGS; i += blockDim.x) {
            geoStringPosXShared[i] = geoStringPosX[i];
            geoStringPosYShared[i] = geoStringPosY[i];
            geoStringMinZShared[i] = geoStringMinZ[i];
            geoStringMaxZShared[i] = geoStringMaxZ[i];
        }
        const float* geoStringPosXPointer = geoStringPosXShared;
        const float* geoStringPosYPointer = geoStringPosYShared;
        const float* geoStringMinZPointer = geoStringMinZShared;
        const float* geoStringMaxZPointer = geoStringMaxZShared;
    #else
        const float* geoStringPosXPointer = geoStringPosX;
        const float* geoStringPosYPointer = geoStringPosY;
        const float* geoStringMinZPointer = geoStringMinZ;
        const float* geoStringMaxZPointer = geoStringMaxZ;
    #endif

    __syncthreads();

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
        bool absorbed = propPhoton(photon, distanceTraveled, rng, scatteringLutPointer, absorbtionLutPointer, absorbtionDeltaTauLutPointer, zOffsetLut);

        // check for collision with DOMs, if collision has happened, the hit will be stored in outputPhotons
        bool collided = checkForCollisionOld(photon, step, distanceTraveled, 
                                  hitIndex, maxHitIndex, outputPhotons, numIndexStringSetPointer, wlenBiasLutPointer, 
                                  geoCellIndex0Pointer, geoCellIndex1Pointer, geoStringInSetPointer, geoLayerNumPointer, geoLayerStartZPointer, geoLayerHeightPointer,
                                  geoStringPosXPointer, geoStringPosYPointer, geoStringMinZPointer, geoStringMaxZPointer);

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

// thread_block_tile::meta_group_rank() requires CUDA 11
#if CUDA_VERSION >= 11000
/**
 * @brief Generates photons for one "step" and simulates propagation through the ice.
 * @param group the group of threads used to process one step (eg one warp)
 * @param step the step to be processed
 * @param sharedPhotonInitials shared memory to store photon initial conditions,
 *                              !! size needs to be same as number of threads in the group !!
 * @param numPhotonsInShared keeps track of the number of photons currently stored in shared memory
 *                           !! needs to live in shared memory itself !!   
 */
__device__ __forceinline__  void propGroup(cg::thread_block_tile<32> group, const I3CLSimStepCuda &step,
                          I3CLPhoton *sharedPhotonInitials, int& numPhotonsInShared,
                          uint32_t *hitIndex, const uint32_t maxHitIndex, I3CLSimPhotonCuda *__restrict__ outputPhotons,
                          const float* wlenLut, const float* zOffsetLut, const float* sharedScatteringLength, 
                          const float* sharedAbsorptionADust, const float* sharedAbsorptionDeltaTau, RngType& rng,
                          const unsigned short *numIndexStringSetPointer,
                          const float* wlenBiasLutPointer,
                          const unsigned short* geoCellIndex0Pointer, 
                          const unsigned short* geoCellIndex1Pointer, 
                          const unsigned char* geoStringInSetPointer, 
                          const unsigned short* geoLayerNumPointer, 
                          const float* geoLayerStartZPointer, 
                          const float* geoLayerHeightPointer,
                          const float* geoStringPosXPointer, 
                          const float* geoStringPosYPointer,  
                          const float* geoStringMinZPointer, 
                          const float* geoStringMaxZPointer)
{
    // calculate step direction
    const float3 stepDir = calculateStepDir(step);

    // variables for managing shared memory
    int photonsLeftInStep = step.numPhotons; // will be 0 or negative if no photons are left

    // local variables for propagating the photon
    int photonId=-1; // threads with a photon id of 0 or bigger contain a valid photon
    I3CLPhoton photon; // this threads current photon

    // generate photon for every thread in the Warp from the step
    if(group.thread_rank() < photonsLeftInStep)
    {
        photon = createPhoton(step, stepDir, wlenLut, rng);
        photonId = 0; // set a valid id
    }
    photonsLeftInStep -= group.size(); // noet: if "photonsLeftInStep" goes negative, it does not matter 

    // make sure shared memory is not in use anymore from previous call
    group.sync();

    // generate photons and store in shared memory
    if(group.thread_rank() < photonsLeftInStep)    
        sharedPhotonInitials[group.thread_rank()] = createPhoton(step, stepDir, wlenLut, rng);
    if(group.thread_rank() == 0) 
    {
        float d = min(group.size(),photonsLeftInStep);
        numPhotonsInShared = d;
        photonsLeftInStep -= d;
    }
    photonsLeftInStep = group.shfl(photonsLeftInStep,0);
    group.sync();
    
    // loop as long as this thread has a valid photon, this is true for all threads while there is a photon left in the "step" or in shared memory
    while(photonId >= 0)
    {
        // propagate photon through the ice
        float distanceTraveled;
        bool absorbed = propPhoton(photon, distanceTraveled, rng, sharedScatteringLength, sharedAbsorptionADust, sharedAbsorptionDeltaTau, zOffsetLut);

        // check for collision with DOMs, if collision has happened, the hit will be stored in outputPhotons
        bool collided = checkForCollisionOld(photon, step, distanceTraveled, 
                                  hitIndex, maxHitIndex, outputPhotons, numIndexStringSetPointer, wlenBiasLutPointer, 
                                  geoCellIndex0Pointer, geoCellIndex1Pointer, geoStringInSetPointer, geoLayerNumPointer, geoLayerStartZPointer, geoLayerHeightPointer,
                                  geoStringPosXPointer, geoStringPosYPointer, geoStringMinZPointer, geoStringMaxZPointer);

        if(collided || absorbed)
        {
            // photon is no longer valid
            photonId = -1;
        }
        else
        {
            // move the photon along its current direction for the distance it was propagated through the ice
            // then scatter to find a new direction vector
            updatePhotonTrack(photon, distanceTraveled);
            scatterPhoton(photon, rng);
        }

        if(numPhotonsInShared > 0 || photonsLeftInStep > 0)
        {
            // there are still photons waiting to be processed as long as this is true, all threads will be in the loop

            if(photonId < 0)
            {
                // try to grab a new photon from shared memory
                photonId = atomicAdd(&numPhotonsInShared,-1)-1;
                if(photonId >= 0)
                    photon = sharedPhotonInitials[photonId];
            }

            // if shared memory is empty, create new photons from the "step" (this branch is taken by all or none of the threads)
            group.sync(); // make sure all threads see the same value of numPhotonsInShared
            if( numPhotonsInShared <= 0 && photonsLeftInStep > 0)
            {
                if(group.thread_rank() < photonsLeftInStep)    
                    sharedPhotonInitials[group.thread_rank()] = createPhoton(step, stepDir, wlenLut, rng);
                    if(group.thread_rank() == 0) 
                    {
                        float d = min(group.size(),photonsLeftInStep);
                        numPhotonsInShared = d;
                        photonsLeftInStep -= d;
                    }
                    photonsLeftInStep = group.shfl(photonsLeftInStep,0);
                    group.sync();

                // if the thread did not have a valid photon, try again to get one now
                if(photonId < 0)
                {
                    // try to grab a new photon from shared memory
                    photonId = atomicAdd(&numPhotonsInShared,-1)-1;
                    if(photonId >= 0)
                        photon = sharedPhotonInitials[photonId];
                }
                group.sync(); // make sure all threads see the same value of numPhotonsInShared for the next iteration
            }
        }
    }
}

/**
 * @brief Main cuda kernel for photon propagation simulation. Called from host with a number of "steps".
 *        Photons for every "step" will be generated, propagated through the ice and stored if they hit the detector.
 *        Does some setup work. Then splits the current work group into thread groups.
 *        Each thread group simulated all photons in one step ( see propGroup() ).
 */
__global__ void propKernelJobqueue(I3CLSimStepCuda* __restrict__ steps, int numSteps, 
                            uint32_t* hitIndex, uint32_t maxHitIndex, I3CLSimPhotonCuda* __restrict__ outputPhotons,
                            const float* wlenLut, const float* zOffsetLut, 
                            uint64_t* __restrict__ rng_x, uint32_t* __restrict__ rng_a, int numPrimes)
{
    #ifdef SHARED_WLEN
        __shared__ float sharedWlenLut[WLEN_LUT_SIZE];
        for (int i = threadIdx.x; i < WLEN_LUT_SIZE; i += blockDim.x) {
            sharedWlenLut[i] = wlenLut[i];
        }
        const float* wlenLutPointer = sharedWlenLut;
    #else
        const float* wlenLutPointer = wlenLut;
    #endif

    #ifdef SHARED_ICE_PROPERTIES
        __shared__ float sharedScatteringLength[171];
        __shared__ float sharedAbsorptionADust[171];
        __shared__ float sharedAbsorptionDeltaTau[171];
        for (int i = threadIdx.x; i < 171; i += blockDim.x) {
            sharedScatteringLength[i] = scatteringLength_b400_LUT[i];
            sharedAbsorptionADust[i] = absorptionLength_aDust400_LUT[i];
            sharedAbsorptionDeltaTau[i] = absorptionLength_deltaTau_LUT[i];
        }
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
        for (int i = threadIdx.x; i < 43; i += blockDim.x) {
            getWavelengthBias_dataShared[i] = getWavelengthBias_data[i];
        }
        const float* wlenBiasLutPointer = getWavelengthBias_dataShared;
    #else
        const float* wlenBiasLutPointer = getWavelengthBias_data;
    #endif

    #ifdef SHARED_NUM_INDEX_STRING_SET
        __shared__ unsigned short geoLayerToOMNumIndexPerStringSetLocal[GEO_geoLayerToOMNumIndexPerStringSet_BUFFER_SIZE];
        for (int i = threadIdx.x; i < GEO_geoLayerToOMNumIndexPerStringSet_BUFFER_SIZE; i += blockDim.x) {
            geoLayerToOMNumIndexPerStringSetLocal[i] = geoLayerToOMNumIndexPerStringSet[i];
        }
        const unsigned short* numIndexStringSetPointer = geoLayerToOMNumIndexPerStringSetLocal;
    #else
        const unsigned short* numIndexStringSetPointer = geoLayerToOMNumIndexPerStringSet;
    #endif

    #ifdef SHARED_COLLISION_GRID_DATA
        __shared__ unsigned short geoCellIndex0shared[GEO_CELL_NUM_X_0 * GEO_CELL_NUM_Y_0];
        __shared__ unsigned short geoCellIndex1shared[GEO_CELL_NUM_X_1 * GEO_CELL_NUM_Y_1];
        for (int i = threadIdx.x; i < GEO_CELL_NUM_X_0 * GEO_CELL_NUM_Y_0; i += blockDim.x) {
            geoCellIndex0shared[i] = geoCellIndex_0[i];
        }
        for (int i = threadIdx.x; i < GEO_CELL_NUM_X_1 * GEO_CELL_NUM_Y_1; i += blockDim.x) {
            geoCellIndex1shared[i] = geoCellIndex_1[i];
        }
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
        for (int i = threadIdx.x; i < NUM_STRINGS; i += blockDim.x) {
            geoStringInSetShared[i] = geoStringInStringSet[i];
        }
        for (int i = threadIdx.x; i < GEO_LAYER_STRINGSET_NUM; i += blockDim.x) {
            geoLayerNumShared[i] = geoLayerNum[i];
            geoLayerStartZShared[i] = geoLayerStartZ[i];
            geoLayerHeightShared[i] = geoLayerHeight[i];
        }
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

    #ifdef SHARED_STRING_POSITIONS
        __shared__ float geoStringPosXShared[NUM_STRINGS];
        __shared__ float geoStringPosYShared[NUM_STRINGS];
        __shared__ float geoStringMinZShared[NUM_STRINGS];
        __shared__ float geoStringMaxZShared[NUM_STRINGS];
        for (int i = threadIdx.x; i < NUM_STRINGS; i += blockDim.x) {
            geoStringPosXShared[i] = geoStringPosX[i];
            geoStringPosYShared[i] = geoStringPosY[i];
            geoStringMinZShared[i] = geoStringMinZ[i];
            geoStringMaxZShared[i] = geoStringMaxZ[i];
        }
        const float* geoStringPosXPointer = geoStringPosXShared;
        const float* geoStringPosYPointer = geoStringPosYShared;
        const float* geoStringMinZPointer = geoStringMinZShared;
        const float* geoStringMaxZPointer = geoStringMaxZShared;
    #else
        const float* geoStringPosXPointer = geoStringPosX;
        const float* geoStringPosYPointer = geoStringPosY;
        const float* geoStringMinZPointer = geoStringMinZ;
        const float* geoStringMaxZPointer = geoStringMaxZ;
    #endif

    __syncthreads();

    // get thread id
    cg::thread_block block = cg::this_thread_block();
    int threadId = blockIdx.x * blockDim.x + threadIdx.x;

    // initialize rng
    RngType rng(rng_x[threadId],rng_a[threadId%numPrimes]);

    // setup shared memory to hold generated photon initial conditions
    __shared__ I3CLPhoton sharedPhotonInitials[NTHREADS_PER_BLOCK];
    __shared__ int numPhotonsInShared[NTHREADS_PER_BLOCK / 32];

    // split up into warps sized groups, each group simulates one step at a time
    cg::thread_block_tile<32> group = cg::tiled_partition<32>(block);
    I3CLPhoton* thisGroupSharedPhotonInitials = &sharedPhotonInitials[0] + (group.size() * group.meta_group_rank());
    const int globalWarpId = blockIdx.x * group.meta_group_size() + group.meta_group_rank();
    const int totalNumWarps = group.meta_group_size() * gridDim.x;

    for(int i = globalWarpId; i < numSteps; i += totalNumWarps)
    {
        const I3CLSimStepCuda step = steps[i];
        propGroup(group, step, thisGroupSharedPhotonInitials, numPhotonsInShared[group.meta_group_rank()], 
                    hitIndex, maxHitIndex, outputPhotons, 
                    wlenLutPointer, zOffsetLut, scatteringLutPointer,
                    absorbtionLutPointer, absorbtionDeltaTauLutPointer, rng,
                    numIndexStringSetPointer, wlenBiasLutPointer, 
                    geoCellIndex0Pointer, geoCellIndex1Pointer, geoStringInSetPointer, geoLayerNumPointer, geoLayerStartZPointer, geoLayerHeightPointer,
                    geoStringPosXPointer, geoStringPosYPointer, geoStringMinZPointer, geoStringMaxZPointer);
    }

    // store rng state
    rng_x[threadId] = rng.getState();
}
#endif // CUDA_VERSION >= 11000
