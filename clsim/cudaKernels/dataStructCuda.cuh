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
    defines all datastructures used in the cuda clsim code
*/

#ifndef DATASTRUCCUDA_CUH
#define DATASTRUCCUDA_CUH

// includes
// ------------------
#include <I3CLSimPhoton.h>
#include <I3CLSimStep.h>
// ------------------

/**
 * @brief Structure holds one simulation "step" to send it from host to device.
 *          Algorithmically this is a work package containing information to create a 
 *          numbermultiple photons. 
 *          Pyhsically, it represents a part of the neutinos path, which is assumed 
 *          to result in  photons with similar properties. 
 */
struct __align__(16) I3CLSimStepCuda {
    float4 posAndTime;           // x,y,z,time             
    float4 dirAndLengthAndBeta;  // theta,phi,length,beta   
    uint32_t numPhotons;         // number of photons created from this step
    uint32_t identifier;
    float weight;                        

    // default constructor, zero everything
    __host__ I3CLSimStepCuda()
        : posAndTime{0, 0, 0, 0},
          dirAndLengthAndBeta{0, 0, 0, 0},
          numPhotons(0),
          identifier(0),
          weight(0)
    {}

    // creates step from the old opencl step structure
    __host__ explicit I3CLSimStepCuda(const I3CLSimStep &i3step)
        : posAndTime({i3step.GetPosX(), i3step.GetPosY(), i3step.GetPosZ(), i3step.GetTime()}),
          dirAndLengthAndBeta({i3step.GetDirTheta(), i3step.GetDirPhi(), i3step.GetLength(), i3step.GetBeta()}),
          numPhotons(i3step.numPhotons),
          identifier(i3step.identifier),
          weight(i3step.weight)
    {}
};

// holds photon when it is written to device memory to be downloaded onto the host
struct __align__(16) I3CLSimPhotonCuda {
    float4 posAndTime;       // x,y,z,time                      
    float2 dir;              // theta,phi                                
    float wavelength;        // photon wavelength                  
    float weight;            
    float groupVelocity;     
    uint32_t identifier;     
    int stringID;          
    int omID;             
    
    __host__ __device__ I3CLSimPhotonCuda() : posAndTime{0, 0, 0, 0} {}

    // load from the original opencl struct
    __host__ I3CLSimPhotonCuda(const I3CLSimPhoton &i3photon)
        : posAndTime{i3photon.GetPosX(), i3photon.GetPosY(), i3photon.GetPosZ(), i3photon.GetTime()},
          dir{i3photon.GetDirTheta(), i3photon.GetDirPhi()},
          wavelength(i3photon.GetWavelength()),
          weight(i3photon.weight),
          identifier(i3photon.GetID()),
          stringID(i3photon.GetStringID()),
          omID(i3photon.GetOMID()),
          groupVelocity(i3photon.GetGroupVelocity())
    {
    }

    // convert to the original opengl struct
    __host__ I3CLSimPhoton getI3CLSimPhoton()
    {
        I3CLSimPhoton i3photon;
        i3photon.SetPosX(posAndTime.x);
        i3photon.SetPosY(posAndTime.y);
        i3photon.SetPosZ(posAndTime.z);
        i3photon.SetTime(posAndTime.w);

        i3photon.SetStartPosX(0);
        i3photon.SetStartPosY(0);
        i3photon.SetStartPosZ(0);
        i3photon.SetStartTime(0);

        i3photon.SetDirTheta(dir.x);
        i3photon.SetDirPhi(dir.y);
        i3photon.SetStartDirTheta(0);
        i3photon.SetStartDirPhi(0);

        i3photon.SetWavelength(wavelength);
        i3photon.SetCherenkovDist(0);
        i3photon.SetNumScatters(0);
        i3photon.SetWeight(weight);
        i3photon.SetID(identifier);
        i3photon.SetStringID(stringID);
        i3photon.SetOMID(omID);
        i3photon.SetGroupVelocity(0);
        i3photon.SetDistInAbsLens(0);

        return i3photon;
    }
};

// holds the initial conditions of a photon in device code
struct I3CLInitialPhoton {
    float3 pos;
    float time;
    float3 dir;
    float wlen;
    float invGroupvel;
    float absLength;
};

// holds photon while it is propagated through the ice
struct I3CLPhoton {
    I3CLPhoton() = default;
    __host__ __device__ explicit I3CLPhoton(const I3CLInitialPhoton &initial)  // creates photon from initial conditions
        : pos(initial.pos),
          time(initial.time),
          dir(initial.dir),
          wlen(initial.wlen),
          invGroupvel(initial.invGroupvel),
          absLength(initial.absLength)
    {}

    float3 pos;
    float time;
    float3 dir;
    float wlen;
    float invGroupvel;
    float absLength;
};

#endif  // DATASTRUCCUDA_CUH
