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

// !! order matters:
#include <propagationKernelSource.cuh>
#include <propagationKernelFunctions.cuh>


 //__device__ global random arrays
 __device__  uint64_t* d_MWC_RNG_x;
 __device__  uint32_t* d_MWC_RNG_a;

 cudaError_t gl_err;
#define CUDA_ERR_CHECK(e) if(cudaError_t(e)!=cudaSuccess) printf("!!! Cuda Error %s in line %d \n", cudaGetErrorString(cudaError_t(e)), __LINE__);
#define CUDA_CHECK_CALL   gl_err = cudaGetLastError(); if(cudaError_t(gl_err)!=cudaSuccess) printf("!!! Cuda Error %s in line %d \n", cudaGetErrorString(cudaError_t(gl_err)), __LINE__-1);
//#define PRINTLC     printf("thread 0 - in line %d \n", __LINE__);  


// remark: ignored tabulate version, removed ifdef TABULATE 
// also removed ifdef DOUBLEPRECISION.
// SAVE_PHOTON_HISTORY  and SAVE_ALL_PHOTONS are not define for now, i.e. commented out these snippets,
// s.t. it corresponds to the default contstructor of I3CLSimStepToPhotonConverterOpenCL

__global__ __launch_bounds__(NTHREADS_PER_BLOCK,4) void
propKernel(uint32_t *hitIndex,   // deviceBuffer_CurrentNumOutputPhotons
      const   uint32_t maxHitIndex, // maxNumOutputPhotons_
#ifndef SAVE_ALL_PHOTONS
const  unsigned short* __restrict__  geoLayerToOMNumIndexPerStringSet,
#endif
  const I3CLSimStepCuda*  __restrict__  inputSteps,      // deviceBuffer_InputSteps
           int nsteps,
           I3CLSimPhotonCuda* __restrict__  outputPhotonsGlobal, // deviceBuffer_OutputPhotons

#ifdef SAVE_PHOTON_HISTORY
           float4 *photonHistory,
#endif
           uint64_t* __restrict__  MWC_RNG_x, uint32_t* __restrict__  MWC_RNG_a);


void photonsToFile(const std::string& filename, I3CLSimPhotonCuda *photons, unsigned int nphotons){
      std::cout<< " writing "<< nphotons << " to file "<< filename<<std::endl;
      std::ofstream outputFile; outputFile.open (filename);
      for (unsigned int i = 0; i <  nphotons; i++)
      {  
                  outputFile <<photons[i].posAndTime.x << "," << photons[i].posAndTime.y << "," << photons[i].numScatters << std::endl;
      }
      outputFile.close();
}


void photonsToFile(const std::string& filename, I3CLSimPhoton *photons, unsigned int nphotons){
      std::cout<< " writing "<< nphotons << " to file "<< filename<<std::endl;
      std::ofstream  outputFile;  outputFile.open (filename);
      for (unsigned int i = 0; i <  nphotons; i++)
      {  
            outputFile <<photons[i].GetPosX() << "," << photons[i].GetPosY() << "," << photons[i].GetNumScatters() << std::endl;
      }
      outputFile.close();
 }

      

// maxNumbWOrkItems from  CL rndm arrays
void init_RDM_CUDA(int maxNumWorkitems, uint64_t* MWC_RNG_x,  uint32_t*  MWC_RNG_a)  
{
	CUDA_ERR_CHECK(cudaMalloc(&d_MWC_RNG_a , maxNumWorkitems* sizeof(uint32_t)));
	CUDA_ERR_CHECK(cudaMalloc(&d_MWC_RNG_x , maxNumWorkitems * sizeof(uint64_t)));
   
	CUDA_ERR_CHECK(cudaMemcpy(d_MWC_RNG_a, MWC_RNG_a,  maxNumWorkitems*sizeof(uint32_t),cudaMemcpyHostToDevice));
      CUDA_ERR_CHECK(cudaMemcpy(d_MWC_RNG_x, MWC_RNG_x, maxNumWorkitems* sizeof(uint64_t),cudaMemcpyHostToDevice));
      
      cudaDeviceSynchronize();
      printf("RNG is set up on CUDA gpu %d. \n", maxNumWorkitems);
}


void launch_CudaPropogate(const I3CLSimStep* __restrict__ in_steps, int nsteps,  
     const uint32_t maxHitIndex, unsigned short *geoLayerToOMNumIndexPerStringSet, int ngeolayer,
        uint64_t* __restrict__  MWC_RNG_x,    uint32_t* __restrict__   MWC_RNG_a, int sizeRNG, float& totalCudaKernelTime, const int nbenchmarks, bool writePhotonsCsv){

       //set up congruental random number generator, reusing host arrays and randomService from I3CLSimStepToPhotonConverterOpenCL setup.
      init_RDM_CUDA( sizeRNG, MWC_RNG_x,  MWC_RNG_a);
      
      printf("nsteps total = %d but dividing into %d launches of max size %d \n", nsteps, 1, nsteps);
      uint32_t h_totalHitIndex =0;
      unsigned short *d_geolayer;
	CUDA_ERR_CHECK(cudaMalloc((void**)&d_geolayer , ngeolayer*sizeof(unsigned short)));
      CUDA_ERR_CHECK(cudaMemcpy(d_geolayer, geoLayerToOMNumIndexPerStringSet, ngeolayer*sizeof(unsigned short),cudaMemcpyHostToDevice));
      
      //these multiple launches correspond to numBuffers..   
      for (int ilaunch= 0 ; ilaunch<1; ++ilaunch )
      {  
            int launchnsteps = ( (1+ilaunch)*nsteps<= nsteps) ?  nsteps : nsteps-ilaunch*nsteps;   
           
            struct I3CLSimStepCuda* h_cudastep = (struct I3CLSimStepCuda*) malloc(launchnsteps*sizeof(struct I3CLSimStepCuda));
            
            for (int i =0; i<launchnsteps; i++){
                  h_cudastep[i] = I3CLSimStep(in_steps[i+ilaunch*nsteps]);      
            } 
                 
            I3CLSimStepCuda * d_cudastep;
            CUDA_ERR_CHECK(cudaMalloc( (void**)&d_cudastep , launchnsteps*sizeof(I3CLSimStepCuda)));
            CUDA_ERR_CHECK(cudaMemcpy(d_cudastep, h_cudastep, launchnsteps*sizeof(I3CLSimStepCuda),cudaMemcpyHostToDevice));
            
            uint32_t *d_hitIndex;
            uint32_t h_hitIndex[1]; h_hitIndex[0]= 0;
            CUDA_ERR_CHECK(cudaMalloc((void**)&d_hitIndex , 1 * sizeof(uint32_t)));
            CUDA_ERR_CHECK(cudaMemcpy(d_hitIndex, h_hitIndex,  1*sizeof(uint32_t),cudaMemcpyHostToDevice));

            I3CLSimPhotonCuda * d_cudaphotons;
            CUDA_ERR_CHECK(cudaMalloc((void**)&d_cudaphotons , maxHitIndex*sizeof(I3CLSimPhotonCuda)));

            int numBlocks =  (launchnsteps+NTHREADS_PER_BLOCK-1)/NTHREADS_PER_BLOCK;
            uint32_t sharedMem =  sizeof(struct I3CLSimStepCuda)*MAX_HITS_PER_SHARED;
 
            const uint32_t maxHitIndexPerThread = maxHitIndex/nsteps;
            printf("remark: avrg over thread maxHitIndex  %u and fixed MAX_HITS_PER_SHARED %d \n",maxHitIndexPerThread,MAX_HITS_PER_SHARED ); 
            printf("launching kernel propKernel<<< %d , %d, %u>>>( .., nsteps=%d)  \n", numBlocks, NTHREADS_PER_BLOCK, sharedMem, launchnsteps);
            propKernel<<<numBlocks, NTHREADS_PER_BLOCK,sharedMem >>>(d_hitIndex, maxHitIndex,  d_geolayer, d_cudastep, launchnsteps, d_cudaphotons, d_MWC_RNG_x, d_MWC_RNG_a);
            cudaDeviceSynchronize(); CUDA_CHECK_CALL

            std::chrono::time_point<std::chrono::system_clock> startKernel = std::chrono::system_clock::now();
            for (int b = 0 ; b< nbenchmarks; ++b){    
                  propKernel<<<numBlocks, NTHREADS_PER_BLOCK,sharedMem>>>(d_hitIndex, maxHitIndex, d_geolayer, d_cudastep, launchnsteps, d_cudaphotons, d_MWC_RNG_x, d_MWC_RNG_a);
            }
              cudaDeviceSynchronize(); 
              std::chrono::time_point<std::chrono::system_clock> endKernel = std::chrono::system_clock::now();
              totalCudaKernelTime =  std::chrono::duration_cast<std::chrono::milliseconds>(endKernel - startKernel).count();
              CUDA_CHECK_CALL
              
            CUDA_ERR_CHECK(cudaMemcpy(h_hitIndex, d_hitIndex, 1*sizeof(uint32_t),cudaMemcpyDeviceToHost));
            uint32_t  numberPhotons = h_hitIndex[0];
            h_totalHitIndex += h_hitIndex[0];
            
            if( numberPhotons> maxHitIndex){
                  printf("Maximum number of photons exceeded, only receiving %" PRIu32
                  " of %" PRIu32 " photons", maxHitIndex, numberPhotons);
                  numberPhotons = maxHitIndex;
            }
     
           // copy (max fo maxHitIndex) photons to host.
            struct I3CLSimPhotonCuda* h_cudaphotons = (struct I3CLSimPhotonCuda*) malloc(numberPhotons*sizeof(struct I3CLSimPhotonCuda));
            CUDA_ERR_CHECK(cudaMemcpy(h_cudaphotons, d_cudaphotons, numberPhotons*sizeof(I3CLSimPhotonCuda),cudaMemcpyDeviceToHost));
            cudaDeviceSynchronize();  

            if(writePhotonsCsv)
            {
                photonsToFile("/home/rhohl/IceCube/offline_production/build/photonsCudaSh.csv",h_cudaphotons, uint32_t(numberPhotons/float(nbenchmarks+1)) );
            }
            
           free(h_cudastep);    
           cudaFree(d_cudaphotons); 
           cudaFree(d_cudastep); 
           cudaFree(d_geolayer); 
      }

      printf( "photon hits = %f from %d steps \n", h_totalHitIndex/float(nbenchmarks+1), nsteps); 
      //check phtoton hits out:
      //for (int i = numberPhotons-10; i< numberPhotons; ++i)printf(" %d photon= %d has hit DOM= %u \n",i, h_cudaphotons[i].stringID , h_cudaphotons[i].omID);       
}


void finalizeCUDA(){
      cudaFree(d_MWC_RNG_a);
      cudaFree(d_MWC_RNG_x);
      cudaDeviceSynchronize();  
}


__global__ void 
propKernel(uint32_t *hitIndex,         // deviceBuffer_CurrentNumOutputPhotons
      const     uint32_t maxHitIndex, // maxNumOutputPhotons_
#ifndef SAVE_ALL_PHOTONS
const   unsigned short* __restrict__ geoLayerToOMNumIndexPerStringSet,
#endif
            const I3CLSimStepCuda* __restrict__ inputSteps,      // deviceBuffer_InputSteps
           int nsteps,
           I3CLSimPhotonCuda* __restrict__ outputPhotonsGlobal, // deviceBuffer_OutputPhotons

#ifdef SAVE_PHOTON_HISTORY
           float4 *photonHistory,
#endif
           uint64_t* __restrict__ MWC_RNG_x, uint32_t* __restrict__ MWC_RNG_a)
           {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  //int global_size = gridDim.x * blockDim.x;

  #ifndef SAVE_ALL_PHOTONS
 
  __shared__    unsigned short   geoLayerToOMNumIndexPerStringSetLocal[GEO_geoLayerToOMNumIndexPerStringSet_BUFFER_SIZE];

for (int ii = threadIdx.x ; ii<GEO_geoLayerToOMNumIndexPerStringSet_BUFFER_SIZE; ii+= blockDim.x){
      geoLayerToOMNumIndexPerStringSetLocal[ii] =geoLayerToOMNumIndexPerStringSet[ii]; 
}  
  __syncthreads();
  
  #endif

 
  //I3CLSimPhotonCuda outputPhotons[MAX_HITS_PER_STEP];
  __shared__ uint32_t localIndexCount, photonsLeftSharedSum, lastThread[2];
  if (threadIdx.x == 0)
  {
      photonsLeftSharedSum = 0;
      localIndexCount = 0;
      lastThread[0] = 0; lastThread[1] =0;
  } 
  extern __shared__  I3CLSimPhotonCuda outputPhotons[MAX_HITS_PER_SHARED ];
    __syncthreads();
 


  if(i >=nsteps) return;
  
 // if(i ==0) printf("CUDA kernel... (thread %u of %u)\n", i, global_size);

  #ifdef SAVE_PHOTON_HISTORY
      float4 currentPhotonHistory[NUM_PHOTONS_IN_HISTORY];
  #endif

  // download MWC RNG state
  uint64_t real_rnd_x = MWC_RNG_x[i];
  uint32_t real_rnd_a = MWC_RNG_a[i];
  uint64_t *rnd_x = &real_rnd_x;
  uint32_t *rnd_a = &real_rnd_a;
 
  
const I3CLSimStepCuda step = inputSteps[i];
  float4 stepDir;
  {
    const float rho = sinf(step.dirAndLengthAndBeta.x); // sin(theta)
    stepDir =
        float4{rho * cosf(step.dirAndLengthAndBeta.y), // rho*cos(phi)
                      rho * sinf(step.dirAndLengthAndBeta.y), // rho*sin(phi)
                      cosf(step.dirAndLengthAndBeta.x),       // cos(phi)
                      ZERO};
  }

 
#define EPSILON 0.00001f
 
  uint32_t photonsLeftToPropagate = step.numPhotons;
  atomicAdd(&photonsLeftSharedSum, 1*(photonsLeftToPropagate<=1) );

      //get thread with highest number of photons
      if( lastThread[1]< photonsLeftToPropagate )
      {
            lastThread[0] = threadIdx.x;
            lastThread[1] = photonsLeftToPropagate;
      }
       
  float abs_lens_left = ZERO;
  float abs_lens_initial = ZERO;

  float4 photonStartPosAndTime;
  float4 photonStartDirAndWlen;
  float4 photonPosAndTime;
  float4 photonDirAndWlen;
  uint32_t photonNumScatters = 0;
  float photonTotalPathLength = ZERO;

#ifdef getTiltZShift_IS_CONSTANT
  int currentPhotonLayer;
#endif

#ifndef FUNCTION_getGroupVelocity_DOES_NOT_DEPEND_ON_LAYER
#error This kernel only works with a constant group velocity (constant w.r.t. layers)
#endif
  float inv_groupvel = ZERO;

  while (photonsLeftToPropagate > 0)
   {
    if (abs_lens_left < EPSILON) {

      // create a new photon
      createPhotonFromTrack(step, stepDir, RNG_ARGS_TO_CALL, photonPosAndTime,
                        photonDirAndWlen);

      // save the start position and time
      photonStartPosAndTime = photonPosAndTime;
      photonStartDirAndWlen = photonDirAndWlen;

      photonNumScatters = 0;
      photonTotalPathLength = ZERO;

#ifdef PRINTF_ENABLED
      printf("   created photon %u at: p=(%f,%f,%f), d=(%f,%f,%f), t=%f, "
             "wlen=%fnm\n", step.numPhotons - photonsLeftToPropagate, photonPosAndTime.x,
             photonPosAndTime.y, photonPosAndTime.z, photonDirAndWlen.x,
             photonDirAndWlen.y, photonDirAndWlen.z, photonPosAndTime.w,
             photonDirAndWlen.w / 1e-9f);
#endif

#ifdef getTiltZShift_IS_CONSTANT
      currentPhotonLayer = min(max(findLayerForGivenZPos(photonPosAndTime.z), 0), MEDIUM_LAYERS - 1);
#endif

      inv_groupvel = 1.f/(getGroupVelocity(0, photonDirAndWlen.w));

      // the photon needs a lifetime. determine distance to next scatter and
      // absorption (this is in units of absorption/scattering lengths)
#ifdef PROPAGATE_FOR_FIXED_NUMBER_OF_ABSORPTION_LENGTHS
      // for table-making, use a fixed number of absorbption lengths
      // (photonics uses a probability of 1e-20, so about 46 absorption lengths)
      abs_lens_initial = PROPAGATE_FOR_FIXED_NUMBER_OF_ABSORPTION_LENGTHS;
#else
      abs_lens_initial = -logf(RNG_CALL_UNIFORM_OC);
#endif
      abs_lens_left = abs_lens_initial;

#ifdef PRINTF_ENABLED
      printf("   - total track length will be %f absorption lengths\n", abs_lens_left);
#endif
    }

    // this block is along the lines of the PPC kernel
    float distancePropagated;
    {
#ifdef getTiltZShift_IS_CONSTANT
#define effective_z (photonPosAndTime.z - getTiltZShift_IS_CONSTANT)
#else
      const float effective_z =  photonPosAndTime.z - getTiltZShift(photonPosAndTime);

      int currentPhotonLayer = min(max(findLayerForGivenZPos(effective_z), 0), MEDIUM_LAYERS - 1);
#endif

      const float photon_dz = photonDirAndWlen.z;

      // add a correction factor to the number of absorption lengths
      // abs_lens_left before the photon is absorbed. This factor will be taken
      // out after this propagation step. Usually the factor is 1 and thus has
      // no effect, but it is used in a direction-dependent way for our model of
      // ice anisotropy.

      const float abs_len_correction_factor = getDirectionalAbsLenCorrFactor(photonDirAndWlen);


      abs_lens_left *= abs_len_correction_factor;

      // the "next" medium boundary (either top or bottom, depending on step
      // direction)
      float mediumBoundary =
          (photon_dz < ZERO) ? (mediumLayerBoundary(currentPhotonLayer))
                             : (mediumLayerBoundary(currentPhotonLayer) +
                                (float)MEDIUM_LAYER_THICKNESS);

      // track this thing to the next scattering point
      float sca_step_left = -logf(RNG_CALL_UNIFORM_OC);
#ifdef PRINTF_ENABLED
      // dbg_printf("   - next scatter in %f scattering lengths\n",
      // sca_step_left);
#endif

 float currentScaLen =  getScatteringLength(currentPhotonLayer, photonDirAndWlen.w);
 float currentAbsLen = getAbsorptionLength(currentPhotonLayer, photonDirAndWlen.w);

  float ais =   (photon_dz * sca_step_left -
           ((mediumBoundary - effective_z))/ currentScaLen) *
          (ONE / (float)MEDIUM_LAYER_THICKNESS);
 float aia =     (photon_dz * abs_lens_left -
           ((mediumBoundary - effective_z))/currentAbsLen) *
          (ONE / (float)MEDIUM_LAYER_THICKNESS);

#ifdef PRINTF_ENABLED
      // dbg_printf("   - ais=%f, aia=%f, j_initial=%i\n", ais, aia,
      // currentPhotonLayer);
#endif

      // propagate through layers
      int j = currentPhotonLayer;
      if (photon_dz < 0) {
        for (; (j > 0) && (ais < ZERO) && (aia < ZERO);
             mediumBoundary -= (float)MEDIUM_LAYER_THICKNESS,
             currentScaLen = getScatteringLength(j, photonDirAndWlen.w),
             currentAbsLen = getAbsorptionLength(j, photonDirAndWlen.w),
             ais += 1.f/(currentScaLen), aia += 1.f/(currentAbsLen))
          --j;
      } else {
        for (; (j < MEDIUM_LAYERS - 1) && (ais > ZERO) && (aia > ZERO);
             mediumBoundary += (float)MEDIUM_LAYER_THICKNESS,
             currentScaLen = getScatteringLength(j, photonDirAndWlen.w),
             currentAbsLen = getAbsorptionLength(j, photonDirAndWlen.w),
             ais -= 1.f/(currentScaLen), aia -= 1.f/(currentAbsLen))
          ++j;
      }

#ifdef PRINTF_ENABLED
      // dbg_printf("   - j_final=%i\n", j);
#endif

      float distanceToAbsorption;
      if ((currentPhotonLayer == j) || ((my_fabs(photon_dz)) < EPSILON)) {
        distancePropagated = sca_step_left * currentScaLen;
        distanceToAbsorption = abs_lens_left * currentAbsLen;
      } else {
        const float recip_photon_dz = 1.f/(photon_dz);
        distancePropagated =
            (ais * ((float)MEDIUM_LAYER_THICKNESS) * currentScaLen +
             mediumBoundary - effective_z) *
            recip_photon_dz;
        distanceToAbsorption =
            (aia * ((float)MEDIUM_LAYER_THICKNESS) * currentAbsLen +
             mediumBoundary - effective_z) *
            recip_photon_dz;
      }
#ifdef getTiltZShift_IS_CONSTANT
      currentPhotonLayer = j;
#endif

#ifdef PRINTF_ENABLED
  //    printf("   - distancePropagated=%f\n", distancePropagated);
#endif

      // get overburden for distance
      if (distanceToAbsorption < distancePropagated) {
        distancePropagated = distanceToAbsorption;
        abs_lens_left = ZERO;
      } else {
        abs_lens_left =
            (distanceToAbsorption - distancePropagated) / currentAbsLen;
      }

      // hoist the correction factor back out of the absorption length
      abs_lens_left = (abs_lens_left)/ abs_len_correction_factor;
    }



#ifndef SAVE_ALL_PHOTONS

    // no photon collission detection in case all photons should be saved
    // the photon is now either being absorbed or scattered.
    // Check for collisions in its way
#ifdef STOP_PHOTONS_ON_DETECTION
#ifdef DEBUG_STORE_GENERATED_PHOTONS
    bool collided;
    if (RNG_CALL_UNIFORM_OC > 0.9) // prescale: 10%
#else                              // DEBUG_STORE_GENERATED_PHOTONS
    bool
#endif                             // DEBUG_STORE_GENERATED_PHOTONS
      
collided =
#endif // STOP_PHOTONS_ON_DETECTION
          checkForCollision(photonPosAndTime, photonDirAndWlen, inv_groupvel,
                            photonTotalPathLength, photonNumScatters,
                            abs_lens_initial - abs_lens_left,
                            photonStartPosAndTime, photonStartDirAndWlen, step,
#ifdef STOP_PHOTONS_ON_DETECTION
                            distancePropagated,
#else  // STOP_PHOTONS_ON_DETECTION
                      distancePropagated,
#endif // STOP_PHOTONS_ON_DETECTION
                        &localIndexCount, maxHitIndex, outputPhotons,
#ifdef SAVE_PHOTON_HISTORY
                            photonHistory, currentPhotonHistory,
#endif // SAVE_PHOTON_HISTORY
                            geoLayerToOMNumIndexPerStringSetLocal);


__syncthreads();
//get thread with highest number of photons
if( lastThread[1]< photonsLeftToPropagate )
{
      lastThread[0] = threadIdx.x;
      lastThread[1] = photonsLeftToPropagate;
}

__syncthreads();

if(localIndexCount >= uint32_t(MAX_HITS_PER_SHARED/2) || photonsLeftSharedSum > 0 ) 
{
    __shared__ uint32_t globStartIndex;

  if(threadIdx.x ==lastThread[0] )
  {
        globStartIndex = atomicAdd(&hitIndex[0], localIndexCount );
  }     
  __syncthreads();
  
  uint32_t sharedIndex = threadIdx.x;

  // ! bug ! todo there is a bug here, because not all threads are active at this point.
    for (; sharedIndex<localIndexCount && globStartIndex+sharedIndex < maxHitIndex ; sharedIndex+= blockDim.x )
  {  
        outputPhotonsGlobal[globStartIndex+sharedIndex] = outputPhotons[sharedIndex]; 
  }
  //reset shared  counter      
  if(threadIdx.x==lastThread[0]){
      localIndexCount = 0;
  } 
  
  __syncthreads();
 
}

__syncthreads();
 

                       
 
#ifdef STOP_PHOTONS_ON_DETECTION
#ifdef DEBUG_STORE_GENERATED_PHOTONS
    collided = true;
#endif // DEBUG_STORE_GENERATED_PHOTONS
    if (collided) {
      // get rid of the photon if we detected it
      abs_lens_left = ZERO;

#ifdef PRINTF_ENABLED
      // dbg_printf("    . colission detected, step limited to
      // thisStepLength=%f!\n",
      //    distancePropagated);
#endif // PRINTF_ENABLED
    }
#endif // STOP_PHOTONS_ON_DETECTION

#endif // not SAVE_ALL_PHOTONS

    // update the track to its next position
    photonPosAndTime.x += photonDirAndWlen.x * distancePropagated;
    photonPosAndTime.y += photonDirAndWlen.y * distancePropagated;
    photonPosAndTime.z += photonDirAndWlen.z * distancePropagated;
    photonPosAndTime.w += inv_groupvel * distancePropagated;
    photonTotalPathLength += distancePropagated;

    // absorb or scatter the photon
    if (abs_lens_left < EPSILON) {
      // photon was absorbed.
      // a new one will be generated at the begin of the loop.
      --photonsLeftToPropagate;
      atomicAdd(&photonsLeftSharedSum, 1*(photonsLeftToPropagate<=1) );
      //get thread with highest number of photons
      if( lastThread[0] == threadIdx.x )
      {
            lastThread[1] = photonsLeftToPropagate;
      }
      


#if defined(SAVE_ALL_PHOTONS) && !defined(TABULATE)
      // save every. single. photon.

      if (RNG_CALL_UNIFORM_CO < SAVE_ALL_PHOTONS_PRESCALE) {
        saveHit(photonPosAndTime, photonDirAndWlen,
                0., // photon has already been propagated to the next position
                inv_groupvel, photonTotalPathLength, photonNumScatters,
                abs_lens_initial, photonStartPosAndTime, photonStartDirAndWlen,
                step,
                0, // string id (not used in this case)
                0, // dom id (not used in this case)
                &localIndexCount, maxHitIndex, outputPhotons
#ifdef SAVE_PHOTON_HISTORY
                ,
                photonHistory, currentPhotonHistory
#endif // SAVE_PHOTON_HISTORY
        );
      }
#endif // SAVE_ALL_PHOTONS 

    } else {
      // photon was NOT absorbed. scatter it and re-start the loop

#ifdef SAVE_PHOTON_HISTORY
      // save the photon scatter point
      currentPhotonHistory[photonNumScatters % NUM_PHOTONS_IN_HISTORY].xyz =
          photonPosAndTime.xyz;
      currentPhotonHistory[photonNumScatters % NUM_PHOTONS_IN_HISTORY].w =
          abs_lens_initial - abs_lens_left;
#endif


// optional direction transformation (for ice anisotropy)
      transformDirectionPreScatter(photonDirAndWlen);

      // choose a scattering angle
      const float cosScatAngle = makeScatteringCosAngle(RNG_ARGS_TO_CALL);
      const float sinScatAngle = sqrt(ONE - sqr(cosScatAngle));

      // change the current direction by that angle
      scatterDirectionByAngle(cosScatAngle, sinScatAngle, photonDirAndWlen,RNG_CALL_UNIFORM_CO);

      // optional direction transformation (for ice anisotropy)
      transformDirectionPostScatter(photonDirAndWlen);
 

      ++photonNumScatters;
  
 
    }
  }

#ifdef PRINTF_ENABLED
  // dbg_printf("Stop kernel... (work item %u of %u)\n", i, global_size);
  // dbg_printf("Kernel finished.\n");
#endif

         // upload MWC RNG state
         MWC_RNG_x[i] = real_rnd_x;
         MWC_RNG_a[i] = real_rnd_a;
        
//copy rest that have not reached MAX_HITS_PER_SHARED

// @todo


}
