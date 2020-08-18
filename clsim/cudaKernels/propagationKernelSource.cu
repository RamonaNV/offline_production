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

#include <propagationKernelSource.cuh>
#include <propagationKernelFunctions.cuh>


#define STATS_TIMERS


cudaError_t gl_err;

#define CUDA_ERR_CHECK(e)              \
    if (cudaError_t(e) != cudaSuccess) \
        printf("!!! Cuda Error %s in line %d \n", cudaGetErrorString(cudaError_t(e)), __LINE__);
#define CUDA_CHECK_CALL                     \
    gl_err = cudaGetLastError();            \
    if (cudaError_t(gl_err) != cudaSuccess) \
        printf("!!! Cuda Error %s in line %d \n", cudaGetErrorString(cudaError_t(gl_err)), __LINE__ - 1);


        
void statsToFile( int* counters, unsigned int n){
    const std::string filename = "/home/rhohl/IceCube/offline_production/build/stats.csv";
    std::cout<< " writing "<< n << " to file "<< filename<<std::endl;
    std::ofstream outputFile; outputFile.open (filename);
    for (unsigned int i = 0; i <  n; i++)
    {  
            outputFile <<counters[i] << std::endl;
    }
    outputFile.close();
}


// remark: ignored tabulate version, removed ifdef TABULATE
// also removed ifdef DOUBLEPRECISION.
// SAVE_PHOTON_HISTORY  and SAVE_ALL_PHOTONS are not define for now, i.e. commented out these snippets,
// s.t. it corresponds to the default contstructor of I3CLSimStepToPhotonConverterOpenCL

__global__ __launch_bounds__(NTHREADS_PER_BLOCK, 4) void propKernel(
    uint32_t* hitIndex,          // deviceBuffer_CurrentNumOutputPhotons
    const uint32_t maxHitIndex,  // maxNumOutputPhotons_
    const unsigned short* __restrict__ geoLayerToOMNumIndexPerStringSet,
    const I3CLSimStepCuda* __restrict__ inputSteps,  // deviceBuffer_InputSteps
    int nsteps,
    I3CLSimPhotonCuda* __restrict__ outputPhotons,  // deviceBuffer_OutputPhotons

#ifdef SAVE_PHOTON_HISTORY
    float4* photonHistory,
#endif
    uint64_t* __restrict__ MWC_RNG_x, uint32_t* __restrict__ MWC_RNG_a
    #ifdef STATS_TIMERS 
    ,clock_t* timers , float* counters, int* perStepCounter
     #endif
    );

// maxNumbWOrkItems from  CL rndm arrays
void init_RDM_CUDA(int maxNumWorkitems, uint64_t* MWC_RNG_x, uint32_t* MWC_RNG_a, uint64_t** d_MWC_RNG_x,
                   uint32_t** d_MWC_RNG_a)
{
    CUDA_ERR_CHECK(cudaMalloc(d_MWC_RNG_a, maxNumWorkitems * sizeof(uint32_t)));
    CUDA_ERR_CHECK(cudaMalloc(d_MWC_RNG_x, maxNumWorkitems * sizeof(uint64_t)));

    CUDA_ERR_CHECK(cudaMemcpy(*d_MWC_RNG_a, MWC_RNG_a, maxNumWorkitems * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_ERR_CHECK(cudaMemcpy(*d_MWC_RNG_x, MWC_RNG_x, maxNumWorkitems * sizeof(uint64_t), cudaMemcpyHostToDevice));

    cudaDeviceSynchronize();
    printf("RNG is set up on CUDA gpu %d. \n", maxNumWorkitems);
}

void launch_CudaPropogate(const I3CLSimStep* __restrict__ in_steps, int nsteps, const uint32_t maxHitIndex,
                          unsigned short* geoLayerToOMNumIndexPerStringSet, int ngeolayer,
                          I3CLSimPhotonSeries& outphotons, uint64_t* __restrict__ MWC_RNG_x,
                          uint32_t* __restrict__ MWC_RNG_a, int sizeRNG, float& totalCudaKernelTime
                     
                        )
{
    #ifdef STATS_TIMERS
      //  nsteps = 1;
    #endif 


    // set up congruental random number generator, reusing host arrays and randomService from
    // I3CLSimStepToPhotonConverterOpenCL setup.
    uint64_t* d_MWC_RNG_x;
    uint32_t* d_MWC_RNG_a;
    init_RDM_CUDA(sizeRNG, MWC_RNG_x, MWC_RNG_a, &d_MWC_RNG_x, &d_MWC_RNG_a);

    printf("nsteps total = %d but dividing into %d launches of max size %d \n", nsteps, 1, nsteps);
    unsigned short* d_geolayer;
    CUDA_ERR_CHECK(cudaMalloc((void**)&d_geolayer, ngeolayer * sizeof(unsigned short)));
    CUDA_ERR_CHECK(cudaMemcpy(d_geolayer, geoLayerToOMNumIndexPerStringSet, ngeolayer * sizeof(unsigned short),
                              cudaMemcpyHostToDevice));
    

    struct I3CLSimStepCuda* h_cudastep = (struct I3CLSimStepCuda*)malloc(nsteps * sizeof(struct I3CLSimStepCuda));

    for (int i = 0; i < nsteps; i++) {
        h_cudastep[i] = I3CLSimStep(in_steps[i]);
        h_cudastep[i].identifier = i;
    }

    #ifdef STATS_TIMERS

   // h_cudastep[0].numPhotons = 1;
    /*
    h_cudastep[0] = I3CLSimStep(in_steps[2649]);
    h_cudastep[0].identifier = 2649;
    for (int i = 0; i < NTHREADS_PER_BLOCK; i++) {
        h_cudastep[i] = I3CLSimStep(in_steps[2649]);
        h_cudastep[i].identifier = 2649;
        h_cudastep[i].numPhotons = 1;
    }
    */
 
    #endif 
    
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

    int numBlocks = (nsteps + NTHREADS_PER_BLOCK - 1) / NTHREADS_PER_BLOCK;
    printf("launching kernel propKernel<<< %d , %d >>>( .., nsteps=%d)  \n", numBlocks, NTHREADS_PER_BLOCK, nsteps);

    #ifdef STATS_TIMERS 
    int ntimers = 5;
    float* d_counters  ;
    clock_t* d_timers  ;
    int * d_perStepCounter;
    const  char *  measurementName[ntimers]; 
    measurementName[1]= "download Photonn and Rng";
    measurementName[2]= "prop Photon";
    measurementName[3]= "check collision";
    measurementName[4]= "while loop";
    measurementName[0] = "whole Kernel";
    
    clock_t timer[numBlocks*2*ntimers];
    float counters[numBlocks*ntimers];
    int perStepCounter[nsteps];

    cudaMalloc((void **)&d_timers, sizeof(clock_t) * ntimers * numBlocks *2);
    cudaMalloc((void **)&d_counters, sizeof(float) * numBlocks*ntimers);
    cudaMalloc((void **)&d_perStepCounter, sizeof(int) * nsteps);

    #endif 


    std::chrono::time_point<std::chrono::system_clock> startKernel = std::chrono::system_clock::now();
    propKernel<<<numBlocks, NTHREADS_PER_BLOCK>>>(d_hitIndex, maxHitIndex, d_geolayer, d_cudastep, nsteps,
                                                  d_cudaphotons, d_MWC_RNG_x, d_MWC_RNG_a
                                                  #ifdef STATS_TIMERS 
                                                  ,d_timers, d_counters,d_perStepCounter
                                                   #endif
                                                );
    CUDA_CHECK_CALL
    CUDA_ERR_CHECK(cudaDeviceSynchronize());
    std::chrono::time_point<std::chrono::system_clock> endKernel = std::chrono::system_clock::now();
    totalCudaKernelTime = std::chrono::duration_cast<std::chrono::milliseconds>(endKernel - startKernel).count();

    CUDA_ERR_CHECK(cudaMemcpy(timer, d_timers, sizeof(clock_t) *ntimers* numBlocks * 2, cudaMemcpyDeviceToHost));
    CUDA_ERR_CHECK(cudaMemcpy(counters, d_counters, sizeof(float) *ntimers* numBlocks, cudaMemcpyDeviceToHost));
    CUDA_ERR_CHECK(cudaMemcpy(perStepCounter,d_perStepCounter , sizeof(int) *nsteps, cudaMemcpyDeviceToHost));

    statsToFile(perStepCounter, nsteps);


       // Compute the difference between the last block end and the first block start.
       clock_t minStart[ntimers];
       clock_t maxEnd[ntimers];
       double avrgTime[ntimers];
       double avrgCounters[ntimers];
 
       for( int m = 0; m< 5; ++m)
       {
           minStart[m] = timer[m* numBlocks + 0];
            maxEnd[m]  = timer[m* numBlocks + numBlocks*ntimers];
            avrgTime[m] = (timer[m* numBlocks + numBlocks*ntimers]-timer[m* numBlocks + 0])/numBlocks; 
            avrgCounters[m] = counters[m* numBlocks + 0]/numBlocks;
            if( m == 2 or m ==3)   avrgTime[m] = double(  timer[m* numBlocks + 0] )/numBlocks; 
       }

       
   
       for (int i = 1; i < numBlocks; i++)
       {
        for( int m = 0; m< 5; ++m)
        {
       
          minStart[m] = timer[m* numBlocks + i] < minStart[m] ? timer[m* numBlocks + i] : minStart[m];
          maxEnd[m] = timer[numBlocks*ntimers+m* numBlocks + i] > maxEnd[m] ? timer[numBlocks*ntimers+m* numBlocks + i] : maxEnd[m];

           avrgTime[m] += double(timer[numBlocks*ntimers+m* numBlocks + i] - timer[m* numBlocks + i] )/numBlocks; 
           avrgCounters[m] += counters[m* numBlocks + i]/numBlocks;
         
           if( m == 2 or m ==3)   avrgTime[m] += double(  timer[m* numBlocks + i] )/numBlocks; 
        }

       }

       printf("clock64 Cycles for one thread, i.e. one step   :  \n" );
       for( int m = 0; m< 5; ++m){
           printf("counted   %f repetition of ' %s '  = %f  and max span %f \n", avrgCounters[m], measurementName[m], avrgTime[m], double(maxEnd[m]-minStart[m]));
       }
           //printf("counted   %f repetition of ' %s '  = %f \n", avrgCounters[m], measurementName[m], avrgTime[m] );
      
       

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
    CUDA_ERR_CHECK(cudaMemcpy(h_cudaphotons, d_cudaphotons, numberPhotons * sizeof(I3CLSimPhotonCuda), cudaMemcpyDeviceToHost));

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
    cudaFree(d_counters);
    cudaFree(d_timers);
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
__device__ __forceinline__ I3CLInitialPhoton createPhoton(const I3CLSimStepCuda &step, float4 stepDir, float* _generateWavelength_0distY, float* _generateWavelength_0distYCumulative, RNG_ARGS)
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
        ph.totalPathLength += distancePropagated;
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

    ++ph.numScatters;
}

__global__ void propKernel(uint32_t* hitIndex,          // deviceBuffer_CurrentNumOutputPhotons
                           const uint32_t maxHitIndex,  // maxNumOutputPhotons_
                           const unsigned short* __restrict__ geoLayerToOMNumIndexPerStringSet,
                           const I3CLSimStepCuda* __restrict__ inputSteps,  // deviceBuffer_InputSteps
                           int nsteps,
                           I3CLSimPhotonCuda* __restrict__ outputPhotons,  // deviceBuffer_OutputPhotons
                           uint64_t* __restrict__ MWC_RNG_x, uint32_t* __restrict__ MWC_RNG_a
                           #ifdef STATS_TIMERS 
                            ,clock_t* timers , float* counters, int* perStepCounter
                            #endif
                        )
{

    #ifdef STATS_TIMERS
        uint64_t start0, end0;
        const int tid = threadIdx.x;
        const int bid = blockIdx.x;

        if(tid == 0){
           for (int m = 0; m<5; ++m) counters[m* gridDim.x + bid] = 0.0;
           for (int m = 0; m<5; ++m) timers[m* gridDim.x + bid] = 0.0;
         
       }
       int m = 0;

        if(tid == 0) start0 = clock64();

      

    #endif
  

#ifndef FUNCTION_getGroupVelocity_DOES_NOT_DEPEND_ON_LAYER
#error This kernel only works with a constant group velocity (constant w.r.t. layers)
#endif

    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

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
    if (i >= nsteps) return;

    
    #ifdef STATS_TIMERS
        uint64_t start, end;
        perStepCounter[i] = 0;
        if(tid == 0) start = clock64();
    #endif

    // download MWC RNG state
    uint64_t real_rnd_x = MWC_RNG_x[i];
    uint32_t real_rnd_a = MWC_RNG_a[i];
    uint64_t* rnd_x = &real_rnd_x;
    uint32_t* rnd_a = &real_rnd_a;

    const I3CLSimStepCuda step = inputSteps[i];

    #ifdef STATS_TIMERS
    if( tid==0 )
    {  
        end = clock64();
        m = 1;
        timers[m* gridDim.x + bid] = start;
        timers[m* gridDim.x + bid + 5*gridDim.x] = end;
        counters[m* gridDim.x + bid] += 1.0;
    }
   #endif


    float4 stepDir;
    {
        const float rho = sinf(step.dirAndLengthAndBeta.x);       // sin(theta)
        stepDir = float4{rho * cosf(step.dirAndLengthAndBeta.y),  // rho*cos(phi)
                         rho * sinf(step.dirAndLengthAndBeta.y),  // rho*sin(phi)
                         cosf(step.dirAndLengthAndBeta.x),        // cos(phi)
                         ZERO};
    }

 
     
    float nphot = step.numPhotons;
    uint32_t photonsLeftToPropagate = step.numPhotons;
    I3CLPhoton photon;
    photon.absLength = 0;
    I3CLInitialPhoton photonInitial;
    
 
    #ifdef STATS_TIMERS
        uint64_t start1, end1;
        if(tid == 0) start1 = clock64();
    #endif


   // #undef STATS_TIMERS
   
    while (photonsLeftToPropagate > 0) {
        
        if (photon.absLength < EPSILON) {
            photonInitial = createPhoton(step, stepDir,_generateWavelength_0distYValuesShared,_generateWavelength_0distYCumulativeValuesShared, RNG_ARGS_TO_CALL);
            photon = I3CLPhoton(photonInitial);
        }

        #ifdef STATS_TIMERS
             perStepCounter[i] += 1;
            if(tid == 0) start = clock64();
        #endif

        // this block is along the lines of the PPC kernel
        float distancePropagated;
        propPhoton(photon, distancePropagated, RNG_ARGS_TO_CALL);

        #ifdef STATS_TIMERS
        if( tid==0 )
        {  
            end = clock64();
            m =2;
            timers[m* gridDim.x + bid] += float(end- start);
            counters[m* gridDim.x + bid] += 1.0;
        }
      #endif


        #ifdef STATS_TIMERS
        if(tid == 0) start = clock64();
        #endif

        bool collided = checkForCollision(photon, photonInitial, step, distancePropagated, 
                                  hitIndex, maxHitIndex, outputPhotons, geoLayerToOMNumIndexPerStringSetLocal, getWavelengthBias_dataShared);

        #ifdef STATS_TIMERS
        if( tid==0 )
        {  
        
            end = clock64();
            m = 3;
            timers[m* gridDim.x + bid] += float(end- start);
            counters[m* gridDim.x + bid] += 1.0;
        }
       #endif

        if (collided) {
            // get rid of the photon if we detected it
            photon.absLength = ZERO;
     //       printf( " photon id = %u of global thread idx %u  and local thread idx %u and block idx %u \n" , photonsLeftToPropagate,i, threadIdx.x, blockIdx.x );
          
        }

        // absorb or scatter the photon
        if (photon.absLength < EPSILON) {
            // photon was absorbed.
            // a new one will be generated at the begin of the loop.
            --photonsLeftToPropagate;
        } else {  // photon was NOT absorbed. scatter it and re-start the loop

            updatePhotonTrack(photon, distancePropagated);
            scatterPhoton(photon, RNG_ARGS_TO_CALL);
        }
    }  // end while
   
    
    
    //#define STATS_TIMERS
    #ifdef STATS_TIMERS
    if( tid==0 )
    {  
       
        end1 = clock64();
        m = 4;
        timers[m* gridDim.x + bid] = start1;
        timers[m* gridDim.x + bid + 5*gridDim.x] = end1;
        counters[m* gridDim.x + bid] += 1.0;

        timers[2* gridDim.x + bid] =timers[2* gridDim.x + bid]/perStepCounter[i]*nphot ;
        timers[3* gridDim.x + bid] =timers[3* gridDim.x + bid]/perStepCounter[i]*nphot ;
    }
  #endif


    // upload MWC RNG state
    MWC_RNG_x[i] = real_rnd_x;
    MWC_RNG_a[i] = real_rnd_a;
    #define STATS_TIMERS
    #ifdef STATS_TIMERS
    if( tid==0 )
    {
          end0 = clock64();
          m = 0;
          timers[m* gridDim.x + bid] = start0;
          timers[m* gridDim.x + bid + 5*gridDim.x] = end0;
          counters[m* gridDim.x + bid] += 1.0;
        
   
    }
    perStepCounter[i] =  int(float(perStepCounter[i])/nphot);

    #endif
 
  
}
