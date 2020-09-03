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

// WARNING: 
// Code corresponds to the default contstructor of I3CLSimStepToPhotonConverterOpenCL
// all extra options have been removed for now
#include <cfenv>

#include <cooperative_groups.h>
#include <cuda/std/atomic>

#include <propagationKernelSource.cuh>
#include <propagationKernelFunctions.cuh>

// namespace alias
namespace cg = cooperative_groups;

cudaError_t gl_err;

#define CUDA_ERR_CHECK(e)              \
    if (cudaError_t(e) != cudaSuccess) \
        printf("!!! Cuda Error %s in line %d \n", cudaGetErrorString(cudaError_t(e)), __LINE__);
#define CUDA_CHECK_CALL                     \
    gl_err = cudaGetLastError();            \
    if (cudaError_t(gl_err) != cudaSuccess) \
        printf("!!! Cuda Error %s in line %d \n", cudaGetErrorString(cudaError_t(gl_err)), __LINE__ - 1);

__global__ __launch_bounds__(NTHREADS_PER_BLOCK, 4) void propKernel(
    uint32_t* hitIndex,          // deviceBuffer_CurrentNumOutputPhotons
    const uint32_t maxHitIndex,  // maxNumOutputPhotons_
    const unsigned short* __restrict__ geoLayerToOMNumIndexPerStringSet,
    const I3CLSimStepCuda* __restrict__ inputSteps,  // deviceBuffer_InputSteps
    int nsteps,
    I3CLSimPhotonCuda* __restrict__ outputPhotons,  // deviceBuffer_OutputPhotons
    uint64_t* __restrict__ MWC_RNG_x, uint32_t* __restrict__ MWC_RNG_a, int numPrimes);

// maxNumbWOrkItems from  CL rndm arrays
void init_RDM_CUDA(int numRngPrimes, int requestedThreads, uint64_t* MWC_RNG_x, uint32_t* MWC_RNG_a, uint64_t** d_MWC_RNG_x,
                   uint32_t** d_MWC_RNG_a)
{
    printf("RNG setup. Available primes: %i. Requested threads %i \n", numRngPrimes, requestedThreads);

    // prime numbers
    CUDA_ERR_CHECK(cudaMalloc(d_MWC_RNG_a, numRngPrimes * sizeof(uint32_t)));
    CUDA_ERR_CHECK(cudaMemcpy(*d_MWC_RNG_a, MWC_RNG_a, numRngPrimes * sizeof(uint32_t), cudaMemcpyHostToDevice));

    // seeds
    std::vector<uint64_t> seeds(requestedThreads);
    uint64_t x = MWC_RNG_x[0];
    int roundMode = std::fegetround();
    std::fesetround(FE_TOWARDZERO);
    for(int i = 0; i < requestedThreads; i++)
    {
        seeds[i]=0;
        while( (seeds[i]==0) | (((uint32_t)(seeds[i]>>32))>=(MWC_RNG_a[i%numRngPrimes]-1)) | (((uint32_t)seeds[i])>=0xfffffffful))
        {
            x = (x & 0xffffffffull) * (MWC_RNG_a[0]) + (x >> 32);
            seeds[i] = x;
        }
    }
    std::fesetround(roundMode);

    CUDA_ERR_CHECK(cudaMalloc(d_MWC_RNG_x, requestedThreads * sizeof(uint64_t)));
    CUDA_ERR_CHECK(cudaMemcpy(*d_MWC_RNG_x, seeds.data(), requestedThreads * sizeof(uint64_t), cudaMemcpyHostToDevice));
}

void launch_CudaPropogate(const I3CLSimStep* __restrict__ in_steps, int nsteps, const uint32_t maxHitIndex,
                          unsigned short* geoLayerToOMNumIndexPerStringSet, int ngeolayer,
                          I3CLSimPhotonSeries& outphotons, uint64_t* __restrict__ MWC_RNG_x,
                          uint32_t* __restrict__ MWC_RNG_a, int sizeRNG, float& totalCudaKernelTime)
{
    int numBlocks = (nsteps*32 + NTHREADS_PER_BLOCK - 1) / NTHREADS_PER_BLOCK; // run 32 threads per step
    int numThreads = numBlocks * NTHREADS_PER_BLOCK;

    // set up congruental random number generator, reusing host arrays and randomService from
    // I3CLSimStepToPhotonConverterOpenCL setup.
    uint64_t* d_MWC_RNG_x;
    uint32_t* d_MWC_RNG_a;
    init_RDM_CUDA(sizeRNG, numThreads, MWC_RNG_x, MWC_RNG_a, &d_MWC_RNG_x, &d_MWC_RNG_a);

    unsigned short* d_geolayer;
    CUDA_ERR_CHECK(cudaMalloc((void**)&d_geolayer, ngeolayer * sizeof(unsigned short)));
    CUDA_ERR_CHECK(cudaMemcpy(d_geolayer, geoLayerToOMNumIndexPerStringSet, ngeolayer * sizeof(unsigned short),
                              cudaMemcpyHostToDevice));

    struct I3CLSimStepCuda* h_cudastep = (struct I3CLSimStepCuda*)malloc(nsteps * sizeof(struct I3CLSimStepCuda));

    for (int i = 0; i < nsteps; i++) {
        h_cudastep[i] = I3CLSimStep(in_steps[i]);
    }

    I3CLSimStepCuda* d_cudastep;
    CUDA_ERR_CHECK(cudaMalloc((void**)&d_cudastep, nsteps * sizeof(I3CLSimStepCuda)));
    CUDA_ERR_CHECK(cudaMemcpy(d_cudastep, h_cudastep, nsteps * sizeof(I3CLSimStepCuda), cudaMemcpyHostToDevice));

    uint32_t* d_hitIndex;
    uint32_t h_hitIndex[1];
    h_hitIndex[0] = 0;
    CUDA_ERR_CHECK(cudaMalloc((void**)&d_hitIndex, 1 * sizeof(uint32_t)));
    CUDA_ERR_CHECK(cudaMemcpy(d_hitIndex, h_hitIndex, 1 * sizeof(uint32_t), cudaMemcpyHostToDevice));

    I3CLSimPhotonCuda* d_cudaphotons;
    CUDA_ERR_CHECK(cudaMalloc((void**)&d_cudaphotons, maxHitIndex * sizeof(I3CLSimPhotonCuda)));

    printf("launching kernel propKernel<<< %d , %d >>>( .., nsteps=%d)  \n", numBlocks, NTHREADS_PER_BLOCK, nsteps);

    std::chrono::time_point<std::chrono::system_clock> startKernel = std::chrono::system_clock::now();
    propKernel<<<numBlocks, NTHREADS_PER_BLOCK>>>(d_hitIndex, maxHitIndex, d_geolayer, d_cudastep, nsteps,
                                                  d_cudaphotons, d_MWC_RNG_x, d_MWC_RNG_a, sizeRNG);

    CUDA_ERR_CHECK(cudaDeviceSynchronize());
    std::chrono::time_point<std::chrono::system_clock> endKernel = std::chrono::system_clock::now();
    totalCudaKernelTime = std::chrono::duration_cast<std::chrono::milliseconds>(endKernel - startKernel).count();

    CUDA_ERR_CHECK(cudaMemcpy(h_hitIndex, d_hitIndex, 1 * sizeof(uint32_t), cudaMemcpyDeviceToHost));
    int numberPhotons = h_hitIndex[0];

    if (numberPhotons > maxHitIndex) {
        printf("Maximum number of photons exceeded, only receiving %" PRIu32 " of %" PRIu32 " photons", maxHitIndex,
               numberPhotons);
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

/**
 * @brief Creates a single photon to be propagated
 * @param step the step to create the photon from
 * @param stepDir step direction to create the photon ( calculated in propGroup() )
 * @param _generateWavelength_0distY data needed for wavelength selection (pass pointer to global or shared data)
 * @param _generateWavelength_0distYCumulative data needed for wavelength selection (pass pointer to global or shared data) 
 * @param RNG_ARGS arguments for the random number generator (use RNG_ARGS_TO_CALL)
 */
__device__ __forceinline__ I3CLInitialPhoton createPhoton(const I3CLSimStepCuda &step, float4 stepDir, const float* _generateWavelength_0distY, const float* _generateWavelength_0distYCumulative, RNG_ARGS)
{
    // create a new photon
    I3CLInitialPhoton ph;
    createPhotonFromTrack(step, stepDir, RNG_ARGS_TO_CALL, ph.posAndTime, ph.dirAndWlen, _generateWavelength_0distY, _generateWavelength_0distYCumulative);
    ph.invGroupvel = 1.f / (getGroupVelocity(0, ph.dirAndWlen.w));

    // set an initial absorption length
    ph.absLength = -logf(RNG_CALL_UNIFORM_OC);
    return ph;
}

/**
 * @brief  propgates a single photon
 * @param ph the photon to propagate
 * @param distancePropagated the distance the photon was propagated during this iteration
 * @param RNG_ARGS arguments for the random number generator (use RNG_ARGS_TO_CALL)
 * @return the propagated distance
 */
__device__ __forceinline__ bool propPhoton(I3CLPhoton& ph, float& distancePropagated, RNG_ARGS)
{ 
    const float effective_z = ph.posAndTime.z - getTiltZShift(ph.posAndTime);
    const int currentPhotonLayer = min(max(findLayerForGivenZPos(effective_z), 0), MEDIUM_LAYERS - 1);
    const float photon_dz = ph.dirAndWlen.z;

    // add a correction factor to the number of absorption lengths
    // abs_lens_left before the photon is absorbed. This factor will be
    // taken out after this propagation step. Usually the factor is 1
    // and thus has no effect, but it is used in a direction-dependent
    // way for our model of ice anisotropy.
    const float abs_len_correction_factor = getDirectionalAbsLenCorrFactor(ph.dirAndWlen);
    ph.absLength *= abs_len_correction_factor;

    // the "next" medium boundary (either top or bottom, depending on
    // step direction)
    float mediumBoundary = (photon_dz < ZERO)
                                ? (mediumLayerBoundary(currentPhotonLayer))
                                : (mediumLayerBoundary(currentPhotonLayer) + (float)MEDIUM_LAYER_THICKNESS);

     // track this thing to the next scattering point
    float scaStepLeft = -logf(RNG_CALL_UNIFORM_OC);

    float currentScaLen = getScatteringLength(currentPhotonLayer, ph.dirAndWlen.w);
    float currentAbsLen = getAbsorptionLength(currentPhotonLayer, ph.dirAndWlen.w);

    float ais = (photon_dz * scaStepLeft - ((mediumBoundary - effective_z)) / currentScaLen) *
                (ONE / (float)MEDIUM_LAYER_THICKNESS);
    float aia = (photon_dz * ph.absLength - ((mediumBoundary - effective_z)) / currentAbsLen) *
                (ONE / (float)MEDIUM_LAYER_THICKNESS);

    
    // propagate through layers
    int j = currentPhotonLayer;
    if (photon_dz < 0) {
        for (; (j > 0) && (ais < ZERO) && (aia < ZERO);
                mediumBoundary -= (float)MEDIUM_LAYER_THICKNESS,
                currentScaLen = getScatteringLength(j, ph.dirAndWlen.w),
                currentAbsLen = getAbsorptionLength(j, ph.dirAndWlen.w), ais += 1.f / (currentScaLen),
                aia += 1.f / (currentAbsLen))
            --j;
    } else {
        for (; (j < MEDIUM_LAYERS - 1) && (ais > ZERO) && (aia > ZERO);
                mediumBoundary += (float)MEDIUM_LAYER_THICKNESS,
                currentScaLen = getScatteringLength(j, ph.dirAndWlen.w),
                currentAbsLen = getAbsorptionLength(j, ph.dirAndWlen.w), ais -= 1.f / (currentScaLen),
                aia -= 1.f / (currentAbsLen))
            ++j;
    }

    float distanceToAbsorption;
    if ((currentPhotonLayer == j) || ((my_fabs(photon_dz)) < EPSILON)) {
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
        ph.absLength = ZERO;
        return true;
    } else {
        ph.absLength = (distanceToAbsorption - distancePropagated) / currentAbsLen;
        
        // hoist the correction factor back out of the absorption length
        ph.absLength = ph.absLength / abs_len_correction_factor;
        return false;
    }

}

/**
 * @brief moves a photon along its track
 * @param ph the photon to move
 * @param distancePropagated the distance the photon was propagated this iteration
 */
__device__ __forceinline__  void updatePhotonTrack(I3CLPhoton& ph, float distancePropagated)
{
        ph.posAndTime.x += ph.dirAndWlen.x * distancePropagated;
        ph.posAndTime.y += ph.dirAndWlen.y * distancePropagated;
        ph.posAndTime.z += ph.dirAndWlen.z * distancePropagated;
        ph.posAndTime.w += ph.invGroupvel * distancePropagated;
}

/**
 * @brief scatters a photon
 * @param ph the photon to scatter
 * @param RNG_ARGS arguments for the random number generator (use RNG_ARGS_TO_CALL) 
 */
__device__ __forceinline__  void scatterPhoton(I3CLPhoton& ph, RNG_ARGS)
{
     // optional direction transformation (for ice anisotropy)
    transformDirectionPreScatter(ph.dirAndWlen);

    // choose a scattering angle
    const float cosScatAngle = makeScatteringCosAngle(RNG_ARGS_TO_CALL);
    const float sinScatAngle = sqrt(ONE - sqr(cosScatAngle));

    // change the current direction by that angle
    scatterDirectionByAngle(cosScatAngle, sinScatAngle, ph.dirAndWlen, RNG_CALL_UNIFORM_CO);

    // optional direction transformation (for ice anisotropy)
    transformDirectionPostScatter(ph.dirAndWlen);
}

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
                          I3CLInitialPhoton *sharedPhotonInitials, int& numPhotonsInShared,
                          uint32_t *hitIndex, const uint32_t maxHitIndex, 
                          I3CLSimPhotonCuda *__restrict__ outputPhotons,
                          const unsigned short* geoLayerToOMNumIndexPerStringSet, const float* _generateWavelength_0distY, 
                          const float* _generateWavelength_0distYCumulative, const float* getWavelengthBias_data,
                          RNG_ARGS)
{
    // calculate step direction
    float4 stepDir;
    const float rho = sinf(step.dirAndLengthAndBeta.x);       // sin(theta)
    stepDir = float4{rho * cosf(step.dirAndLengthAndBeta.y),  // rho*cos(phi)
                        rho * sinf(step.dirAndLengthAndBeta.y),  // rho*sin(phi)
                        cosf(step.dirAndLengthAndBeta.x),        // cos(phi)
                        ZERO};
    

    // variables for managing shared memory
    int photonsLeftInStep = step.numPhotons; // will be 0 or negative if no photons are left

    // local variables for propagating the photon
    int photonId=-1; // threads with a photon id of 0 or bigger contain a valid photon
    I3CLPhoton photon; // this threads current photon

    // generate photon for every thread in the Warp from the step
    if(group.thread_rank() < photonsLeftInStep)
    {
        I3CLInitialPhoton photonInitial = createPhoton(step, stepDir, _generateWavelength_0distY, _generateWavelength_0distYCumulative, RNG_ARGS_TO_CALL);
        photon = I3CLPhoton(photonInitial);
        photonId = 0; // set a valid id
    }
    photonsLeftInStep -= group.size(); // noet: if "photonsLeftInStep" goes negative, it does not matter 

    // make sure shared memory is not in use anymore from previous call
    group.sync();

    // generate photons and store in shared memory
    if(group.thread_rank() < photonsLeftInStep)    
        sharedPhotonInitials[group.thread_rank()] = createPhoton(step, stepDir, _generateWavelength_0distY, _generateWavelength_0distYCumulative, RNG_ARGS_TO_CALL);
    if(group.thread_rank() == 0) 
        numPhotonsInShared = min(group.size(),photonsLeftInStep);
    group.sync();
    photonsLeftInStep -= numPhotonsInShared;
    group.sync(); // make sure numPhotonsInShared is not changed before photonsLeftInStep was updated (maybe move photonsLeftInStep to shared?)
    
    // loop as long as this thread has a valid photon, this is true for all threads while there is is photon left in the "step" or in shared memory
    while(photonId >= 0)
    {
        // propagate photon through the ice
        float distancePropagated;
        bool absorbed = propPhoton(photon, distancePropagated, RNG_ARGS_TO_CALL);

        // check if a collision with the sensor did occur
        bool collided = checkForCollision(photon, step, distancePropagated, 
                                  hitIndex, maxHitIndex, outputPhotons, geoLayerToOMNumIndexPerStringSet, getWavelengthBias_data);

        if(collided)
        {
            // TODO: store (currently it is stored inside the collision check which is very confusing)
            // photon is no longer valid
            photonId = -1;
        }
        else if(absorbed)
        {
            // photon is no longer valid
            photonId = -1;
        }
        else
        {
            // move the photon, scatter and repeat
            updatePhotonTrack(photon, distancePropagated);
            scatterPhoton(photon, RNG_ARGS_TO_CALL);
        }

        if(numPhotonsInShared > 0 || photonsLeftInStep > 0)
        {
            // there are still photons waiting to be processed as long as this is true, all threads will be in the loop

            if(photonId < 0)
            {
                // try to grab a new photon from shared memory
                photonId = atomicAdd(&numPhotonsInShared,-1)-1;
                if(photonId >= 0)
                {
                    I3CLInitialPhoton photonInitial = sharedPhotonInitials[photonId];
                    photon = I3CLPhoton(photonInitial);
                }
            }

            // if shared memory is empty, create new photons from the "step" (this branch is taken by all or none of the threads)
            group.sync(); // make sure all threads see the same value of numPhotonsInShared
            if( numPhotonsInShared <= 0 && photonsLeftInStep > 0)
            {
                if(group.thread_rank() < photonsLeftInStep)    
                    sharedPhotonInitials[group.thread_rank()] = createPhoton(step, stepDir, _generateWavelength_0distY, _generateWavelength_0distYCumulative, RNG_ARGS_TO_CALL);
                if(group.thread_rank() == 0) 
                    numPhotonsInShared = min(group.size(),photonsLeftInStep);
                group.sync();
                photonsLeftInStep -= numPhotonsInShared;
                group.sync(); // make sure numPhotonsInShared is not changed before photonsLeftInStep was updated (maybe move photonsLeftInStep to shared?)

                // if the thread did not have a valid photon, try again to get one now
                if(photonId < 0)
                {
                    // try to grab a new photon from shared memory
                    photonId = atomicAdd(&numPhotonsInShared,-1)-1;
                    if(photonId >= 0)
                    {
                        I3CLInitialPhoton photonInitial = sharedPhotonInitials[photonId];
                        photon = I3CLPhoton(photonInitial);
                    }    
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
__global__ void propKernel(uint32_t *hitIndex, const uint32_t maxHitIndex,
                           const unsigned short *__restrict__ geoLayerToOMNumIndexPerStringSet,
                           const I3CLSimStepCuda *__restrict__ inputSteps, int nsteps,
                           I3CLSimPhotonCuda *__restrict__ outputPhotons,
                           uint64_t *__restrict__ MWC_RNG_x, uint32_t *__restrict__ MWC_RNG_a, int numPrimes)
{
    cg::thread_block block = cg::this_thread_block();
    int threadId = blockIdx.x * blockDim.x + threadIdx.x;

    // cache some data in shared memory
    __shared__ unsigned short geoLayerToOMNumIndexPerStringSetLocal[GEO_geoLayerToOMNumIndexPerStringSet_BUFFER_SIZE];
    __shared__ float _generateWavelength_0distYValuesShared[_generateWavelength_0NUM_DIST_ENTRIES];
    __shared__ float _generateWavelength_0distYCumulativeValuesShared[_generateWavelength_0NUM_DIST_ENTRIES];
    __shared__ float getWavelengthBias_dataShared[_generateWavelength_0NUM_DIST_ENTRIES];

    for (int ii = threadIdx.x; ii < GEO_geoLayerToOMNumIndexPerStringSet_BUFFER_SIZE; ii += blockDim.x) {
        geoLayerToOMNumIndexPerStringSetLocal[ii] = geoLayerToOMNumIndexPerStringSet[ii];
    }
    for (int ii = threadIdx.x; ii < _generateWavelength_0NUM_DIST_ENTRIES; ii += blockDim.x) {
        _generateWavelength_0distYValuesShared[ii] = _generateWavelength_0distYValues[ii];
        _generateWavelength_0distYCumulativeValuesShared[ii] = _generateWavelength_0distYCumulativeValues[ii];
        getWavelengthBias_dataShared[ii] = getWavelengthBias_data[ii];
    }
    __syncthreads();

    // download MWC RNG state
    uint64_t real_rnd_x = MWC_RNG_x[threadId];
    uint32_t real_rnd_a = MWC_RNG_a[threadId%numPrimes];
    uint64_t *rnd_x = &real_rnd_x;
    uint32_t *rnd_a = &real_rnd_a;

    // setup shared memory to hold generated photon initial conditions
    __shared__ I3CLInitialPhoton sharedPhotonInitials[NTHREADS_PER_BLOCK];
    __shared__ int numPhotonsInShared[NTHREADS_PER_BLOCK / 32];

    // split up into warps sized groups, each group simulates one step at a time
    cg::thread_block_tile<32> group = cg::tiled_partition<32>(block);
    I3CLInitialPhoton* thisGroupSharedPhotonInitials = &sharedPhotonInitials[0] + (group.size() * group.meta_group_rank());
    const int globalWarpId = blockIdx.x * group.meta_group_size() + group.meta_group_rank();
    const int totalNumWarps = group.meta_group_size() * gridDim.x;

    for(int i = globalWarpId; i < nsteps; i += totalNumWarps)
    {
        const I3CLSimStepCuda step = inputSteps[i];
        propGroup(group, step, thisGroupSharedPhotonInitials, numPhotonsInShared[group.meta_group_rank()], 
                    hitIndex, maxHitIndex, outputPhotons, 
                    geoLayerToOMNumIndexPerStringSetLocal, _generateWavelength_0distYValuesShared,
                    _generateWavelength_0distYCumulativeValuesShared, getWavelengthBias_dataShared, 
                    RNG_ARGS_TO_CALL);
    }

    // upload MWC RNG state
    MWC_RNG_x[threadId] = real_rnd_x;
}