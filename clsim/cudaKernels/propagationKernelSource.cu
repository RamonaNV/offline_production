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

#include <propagationKernelSource.cuh>
#include <propagationKernelFunctions.cuh>
 
#include <chrono>
#include <string>
#include <vector>
 
#include <rng.cuh>
 
#include <optixSrc/common.h>
#include <optixSrc/params.hpp>
#include <optixSrc/rtxFunctions.hpp>
#include <optixSrc/io_utils.hpp>
#include <optix.h>

RTXDataHolder *rtx_dataholder; 

void launch_CudaPropogate(const I3CLSimStep* __restrict__ in_steps, uint32_t nsteps, const uint32_t maxHitIndex,
                          unsigned short* geoLayerToOMNumIndexPerStringSet, int ngeolayer,
                          I3CLSimPhotonSeries& outphotons, uint64_t* __restrict__ MWC_RNG_x,
                          uint32_t* __restrict__ MWC_RNG_a, int sizeRNG, float& totalCudaKernelTime)
{

  // ---------------------------- set up Optix -------------------------------------

 
  std::string domfile =   std::string(DATA_DIR) + "doms.obj"; //pillDOMs doms
  std::string stringfile =   std::string(DATA_DIR) + "strings.obj";
  std::cout << "DOM file = " << domfile << std::endl;
  std::cout << "String file = " << stringfile << std::endl;
  std::string ptx_filename = PTX_DIR "optixKernels.ptx";

  std::vector<std::string> obj_files;
  //order matters, addd doms, then strings
  obj_files.push_back(domfile);
  bool simulate_cables = true   ; // todo add to parameter file
  if(simulate_cables)    obj_files.push_back(stringfile);

  cudaStream_t stream;
  CUDA_CHECK(cudaStreamCreate(&stream));

  rtx_dataholder = new RTXDataHolder();
  if(OPTIX_VERBOSE) std::cout << "Initializing Context \n";
  rtx_dataholder->initContext();
  if(OPTIX_VERBOSE)   std::cout << "Reading PTX file and creating modules \n";
  rtx_dataholder->createModule(ptx_filename);
  if(OPTIX_VERBOSE)  std::cout << "Creating Optix Program Groups \n";
  rtx_dataholder->createProgramGroups();
  if(OPTIX_VERBOSE)  std::cout << "Linking Pipeline \n";
  rtx_dataholder->linkPipeline();
  if(OPTIX_VERBOSE)  std::cout << "Building Shader Binding Table (SBT) \n";
  rtx_dataholder->buildSBT();


  if(OPTIX_VERBOSE)  std::cout << "Building Acceleration Structure \n";
  bool setGAS  = rtx_dataholder->buildAccelerationStructure(obj_files);

  // ---------------------------- preapare device pointers -------------------------------------

  // set up congruental random number generator, reusing host arrays and randomService from
  // I3CLSimStepToPhotonConverterOpenCL setup.
  uint64_t* d_MWC_RNG_x;
  uint32_t* d_MWC_RNG_a;
  initMWCRng(sizeRNG, MWC_RNG_x, MWC_RNG_a, &d_MWC_RNG_x, &d_MWC_RNG_a);

  struct I3CUDASimStep* h_cudastep = (struct I3CUDASimStep*)malloc(nsteps * sizeof(struct I3CUDASimStep));

  unsigned long countPhotons = 0;
  for (uint32_t i = 0; i < nsteps; i++) {
    h_cudastep[i] = I3CUDASimStep(in_steps[i]);
    countPhotons += in_steps[i].numPhotons;
  }

  I3CUDASimStep* d_cudastep;
  CUDA_ERR_CHECK(cudaMalloc((void**)&d_cudastep, nsteps * sizeof(I3CUDASimStep)));
  CUDA_ERR_CHECK(cudaMemcpy(d_cudastep, h_cudastep, nsteps * sizeof(I3CUDASimStep), cudaMemcpyHostToDevice));

  uint32_t *d_hitIndex;
  CUDA_CHECK(cudaMalloc((void **)&d_hitIndex, 1* sizeof(uint32_t)));
  CUDA_CHECK(cudaMemset(d_hitIndex, 0,1*sizeof(uint32_t)));

  uint32_t *d_cableHits;
  CUDA_CHECK(cudaMalloc((void **)&d_cableHits, 1* sizeof(uint32_t)));
  CUDA_CHECK(cudaMemset(d_cableHits, 0,1*sizeof(uint32_t)));


  I3CUDASimPhoton* d_photons;
  CUDA_ERR_CHECK(cudaMalloc((void**)&d_photons, maxHitIndex * sizeof(I3CUDASimPhoton)));

  auto wlenLut = generateWavelengthLut();
  float* d_wlenLut;
  CUDA_ERR_CHECK(cudaMalloc((void**)&d_wlenLut, WLEN_LUT_SIZE * sizeof(float)));
  CUDA_ERR_CHECK(cudaMemcpy(d_wlenLut, wlenLut.data(), WLEN_LUT_SIZE * sizeof(float), cudaMemcpyHostToDevice));

  auto zOffsetLut = generateZOffsetLut();
  float* d_zOffsetLut;
  CUDA_ERR_CHECK(cudaMalloc((void**)&d_zOffsetLut, zOffsetLut.size() * sizeof(float)));
  CUDA_ERR_CHECK(cudaMemcpy(d_zOffsetLut, zOffsetLut.data(), zOffsetLut.size() * sizeof(float), cudaMemcpyHostToDevice));

  // ---------------------------- add device pointers to Optix Parameters -------------------------------------
  Params params;
  params.handle = rtx_dataholder->gas_handle;
  params.nsteps = nsteps;
  params.maxHitIndex =  maxHitIndex;
  params.cableHits =  d_cableHits;
  params.hitIndex = d_hitIndex;
  params.steps = d_cudastep;
  params.MWC_RNG_x = d_MWC_RNG_x;
  params.MWC_RNG_a = d_MWC_RNG_a;
  params.hitPhotons = d_photons;
  params.wlenLut = d_wlenLut;
  params.zOffsetLut = d_zOffsetLut;

  Params *d_param;
  CUDA_CHECK(cudaMalloc((void **)&d_param, sizeof(Params)));
  CUDA_CHECK(cudaMemcpy(d_param, &params, sizeof(params), cudaMemcpyHostToDevice));

  const OptixShaderBindingTable &sbt = rtx_dataholder->sbt;


  // ---------------------------- launching Optix -------------------------------------
  std::cout << "Launching OptiX for "<< nsteps <<" steps \n"; 

  //prerun
  OPTIX_CHECK(optixLaunch(rtx_dataholder->pipeline, stream,
    reinterpret_cast<CUdeviceptr>(d_param),
    sizeof(Params), &sbt, nsteps, 1, 1));
  CUDA_CHECK(cudaDeviceSynchronize());

  const int nruns = 5;
  auto beg = std::chrono::steady_clock::now();
  for (int i = 0 ; i<nruns; i++)
  {
    OPTIX_CHECK(optixLaunch(rtx_dataholder->pipeline, stream,
                            reinterpret_cast<CUdeviceptr>(d_param),
                            sizeof(Params), &sbt, nsteps, 1, 1));
  }
  
  CUDA_CHECK(cudaDeviceSynchronize());
  auto end = std::chrono::steady_clock::now();
  double timer = std::chrono::duration_cast<std::chrono::duration<double>>((end-beg)).count();


  std::cout << " --------------------------------- \n OptiX kernel time, averaged over "<< nruns << " runs = "
    << timer*1000/nruns << " [ms] \n --------------------------------- "<< std::endl;

    std::cout << " --------------------------------- \n OptiX kernel time, averaged over "<< nruns << " runs = "
    << timer/countPhotons*1000000000/nruns << " [ns/photon] \n --------------------------------- "<< std::endl;



  // ---------------------------- copy back to host  ------------------------------------- 
  uint32_t numberHits = 0;
  CUDA_CHECK(cudaMemcpy(&numberHits, d_hitIndex , 1*sizeof(uint32_t), cudaMemcpyDeviceToHost));


  uint32_t cableHits = 0;
  CUDA_CHECK(cudaMemcpy(&cableHits, d_cableHits , 1*sizeof(uint32_t), cudaMemcpyDeviceToHost));


  numberHits /= (nruns+1.0);
  cableHits  /= (nruns+1.0);
  std::cout<<" % of hits = "<< 100.0*double(numberHits)/countPhotons<< ", number of DOM hits = "<<  numberHits <<", number of cable hits = "<<  cableHits <<", number of tot hits = "<< numberHits+cableHits <<", number photons = "<<countPhotons << std::endl;  

  if (numberHits > maxHitIndex) {
      std::cout<< "Maximum number of photons exceeded, only receiving "  <<maxHitIndex <<" out of " << numberHits<<std::endl; 
  }

   // copy (max fo maxHitIndex) photons to host
    struct I3CUDASimPhoton* h_photons =
        (struct I3CUDASimPhoton*)malloc(numberHits * sizeof(struct I3CUDASimPhoton));
    CUDA_ERR_CHECK(cudaMemcpy(h_photons, d_photons, numberHits * sizeof(I3CUDASimPhoton), cudaMemcpyDeviceToHost));
     
    outphotons.resize(numberHits);

    //for testing
    float3* hitPositions = (float3*) std::malloc(numberHits*sizeof(float3));
    if(TEST_FEASABLE_HIT_POS || WRITE_OUTPHOTONS){
       for (int i = 0; i < numberHits; i++) {
          hitPositions[i] = float3{h_photons[i].posAndTime.x, h_photons[i].posAndTime.y,h_photons[i].posAndTime.z};
          outphotons[i] = h_photons[i].getI3CLSimPhoton();
      }
    }

// ----------------------------  verification / to check if results are feasible  -------------------------------------

if(WRITE_OUTPHOTONS) photonsToFile("./outphotons_optix.csv", h_photons, numberHits) ;
if(TEST_FEASABLE_HIT_POS) verifyHitsLocations(hitPositions, int(numberHits) );
 
// ----------------------------  clean up -------------------------------------
if(OPTIX_VERBOSE)  std::cout << "Cleaning up ...  \n";
  free(h_cudastep);
  free(h_photons);
  CUDA_CHECK(cudaFree(d_photons));
  CUDA_CHECK( cudaFree(d_cudastep));
  CUDA_CHECK( cudaFree(d_MWC_RNG_a));
  CUDA_CHECK( cudaFree(d_MWC_RNG_x));
  CUDA_CHECK( cudaFree(d_wlenLut));
  CUDA_CHECK( cudaFree(d_zOffsetLut));
  CUDA_CHECK(cudaFree(d_param));
  delete rtx_dataholder;
  std::cout <<" ------------------ optiX - done --------------------- \n \n";
}

