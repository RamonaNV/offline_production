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


#include <optix.h> 
#include<utils.cuh>

#include <params.hpp>
#include <rng.cuh>
#include <propagationKernelFunctions.cuh>



 extern "C" static __constant__ Params params;
  
 /*
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
__device__ __forceinline__ bool propPhoton(I3CUDAPhoton& ph, float& distanceTraveled, RngType& rng, const float* scatteringLength, const float* absorptionDust, const float* absorptionDeltaTau, const float* zOffsetLut);

 
__device__ __forceinline__ bool propPhoton(I3CUDAPhoton& ph, float& distancePropagated, RngType& rng, const float* scatteringLength, const float* absorptionDust, const float* absorptionTauDelta, const float* zOffsetLut)
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
 
 extern "C" __global__ void __raygen__prog() {
  const uint32_t launch_index_x = optixGetLaunchIndex().x;
   
    // download MWC RNG state
    RngType rng(params.MWC_RNG_x[launch_index_x], params.MWC_RNG_a[launch_index_x]);

    const float* wlenLutPointer = params.wlenLut;
    //const float* wlenBiasLutPointer = getWavelengthBias_data;
    const float* scatteringLutPointer = scatteringLength_b400_LUT;
    const float* absorbtionLutPointer = absorptionLength_aDust400_LUT;
    const float* absorbtionDeltaTauLutPointer = absorptionLength_deltaTau_LUT;

    const I3CUDASimStep step = params.steps[launch_index_x];
    const float3 stepDir = detail::calculateStepDir(step);

    // variables to store data about current photon
    uint32_t photonsLeftToPropagate = step.numPhotons;
    I3CUDAPhoton photon;
    photon.absLength = 0.0f;

    // loop until all photons are done
    while (photonsLeftToPropagate > 0) {

        if (photon.absLength < EPSILON) {
                  photon = detail::createPhoton(step, stepDir, wlenLutPointer, rng);
        }

        // this block is along the lines of the PPC kernel
        float distancePropagated;
        bool absorbed = propPhoton(photon, distancePropagated, rng, scatteringLutPointer, absorbtionLutPointer, absorbtionDeltaTauLutPointer, params.zOffsetLut);


        // ----------------- check for collision part, i.e. trace ray ---------------------
           
        // create a ray
        float3 ray_origin;
        float3 ray_direction;
      
        ray_origin.x = photon.pos.x;
        ray_origin.y = photon.pos.y;
        ray_origin.z = photon.pos.z;

        ray_direction.x = photon.dir.x;
        ray_direction.y = photon.dir.y;
        ray_direction.z = photon.dir.z;
       
        float tmin = 0.0f;
        float tmax = distancePropagated;
      
        float ray_time = 0.0f;
        OptixVisibilityMask visibilityMask = 255;
        const unsigned int rayFlags = OPTIX_RAY_FLAG_DISABLE_ANYHIT;
        const  unsigned int SBToffset = 0;
        const  unsigned int SBTstride = 0;
        const  unsigned int missSBTIndex = 0;
      
        // Extract Payload as unsigned int
        unsigned int hit_type = MISS;
      
        unsigned int wlen = float_as_uint(photon.wlen);
        unsigned int groupvel = float_as_uint(photon.invGroupvel);
        unsigned int time = float_as_uint(photon.wlen); //todo

        optixTrace(params.handle, ray_origin, ray_direction, tmin, tmax, ray_time,
                    visibilityMask, rayFlags, SBToffset, SBTstride, missSBTIndex,
                    hit_type,  wlen,groupvel,time );          
    
      if(hit_type != MISS )
      {
        photon.absLength = 0.0;
      }  

       // absorb or scatter the photon
        if ( photon.absLength < EPSILON) {
            // photon was absorbed.
            // a new one will be generated at the begin of the loop.
            --photonsLeftToPropagate;
        } else {  // photon was NOT absorbed. scatter it and re-start the loop

            detail::updatePhotonTrack(photon, distancePropagated);
            detail::scatterPhoton(photon, rng);
        }
    }  // end while

    // upload MWC RNG state
   params.MWC_RNG_x[launch_index_x] = rng.getState();
 }
 
 
 extern "C" __global__ void __closesthit__prog_DOMs() {

        optixSetPayload_0(DOM);

         const uint32_t launch_index_x = optixGetLaunchIndex().x;
        unsigned int tri_id = optixGetPrimitiveIndex();
        // We defined out geometry as a triangle geometry. In this case the
        // We add the t value of the intersection
        float ray_tmax =  optixGetRayTmax();

        float3 ray_dir = optixGetWorldRayDirection();
        float3 ray_origin = optixGetWorldRayOrigin();
        float3 hit_point = ray_origin + ray_tmax * ray_dir;

        // printf("Hit at triangle = %u, hit point = ( %f,  %f, %f)\n", tri_id, hit_point.x, hit_point.y, hit_point.z);
        const uint32_t myIndex =  atomicAdd(params.hitIndex, 1);

        // this seems to take a lot of time, mostly filling the photon 
        if(myIndex < params.maxHitIndex )
        {       
            I3CUDASimPhoton outphoton;
            const float distance = sqrtf( (ray_origin.x-hit_point.x)*(ray_origin.x-hit_point.x) + (ray_origin.y-hit_point.y)*(ray_origin.y-hit_point.y) + (ray_origin.z-hit_point.z)*(ray_origin.z- hit_point.z) ); 
            outphoton.posAndTime = float4{hit_point.x,hit_point.y,hit_point.z, uint_as_float( optixGetPayload_3() ) + distance * uint_as_float( optixGetPayload_2() ) };
            outphoton.dir = detail::sphDirFromCar(ray_dir); 
            outphoton.wavelength = uint_as_float( optixGetPayload_1() );  
            outphoton.omID = ushort(tri_id); //todo
            outphoton.weight =  params.steps[launch_index_x].weight/ getWavelengthBias(uint_as_float( optixGetPayload_1() ), getWavelengthBias_data); 
            outphoton.identifier = params.steps[launch_index_x].identifier ;  
            outphoton.groupVelocity = 1.f /  uint_as_float( optixGetPayload_2() );  

            //unnecessary or not needed according to icecube guys
            /*  outphoton.stringID = short(tri_id);
            outphoton.cherenkovDist = photon.totalPathLength + distance;
            outphoton.numScatters = photon.numScatters;
            outphoton.startPosAndTime = photonInitial.posAndTime;  
            outphoton.startDir = sphDirFromCar(photonInitial.dirAndWlen);  
            outphoton.distInAbsLens = photonInitial.absLength - photon.absLength;   */

            params.hitPhotons[myIndex] = outphoton;
        } 

       // printf( " hit ");

}  

/*extern "C" __global__ void __closesthit__prog_cables() {    
    atomicAdd(params.cableHits, 1);
    optixSetPayload_0(CABLE);
}  
 */