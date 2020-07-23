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

#ifndef DATASTRUCCUDA_CUH
#define DATASTRUCCUDA_CUH

#include <I3CLSimPhoton.h>
#include <I3CLSimStep.h>
struct I3CLSimStepCuda {
    float4 posAndTime;           // x,y,z,time                      // 4x  4byte float
    float4 dirAndLengthAndBeta;  // theta,phi,length,beta    // 4x 4byte float
    uint32_t numPhotons;         //    4byte unsigned
    float weight;                //    4byte float
    uint32_t identifier;         //    4byte unsigned
    unsigned char sourceType;    //     1 byte unsigned
                                 // total: 45 bytes

    __host__ I3CLSimStepCuda()
        : posAndTime{0, 0, 0, 0},
          dirAndLengthAndBeta{0, 0, 0, 0},
          numPhotons(0),
          weight(0),
          identifier(0),
          sourceType(0){};

    __host__ I3CLSimStepCuda &operator=(const I3CLSimStep &i3step)
    {
        posAndTime = float4{i3step.GetPosX(), i3step.GetPosY(), i3step.GetPosZ(), i3step.GetTime()};
        dirAndLengthAndBeta = float4{i3step.GetDirTheta(), i3step.GetDirPhi(), i3step.GetLength(), i3step.GetBeta()};
        numPhotons = i3step.numPhotons;
        weight = i3step.weight;
        identifier = i3step.identifier;
        sourceType = i3step.sourceType;

        return *this;
    }
};

struct I3CLSimPhotonCuda {
    float4 posAndTime;       // x,y,z,time                      // 4x 32bit float
    float2 dir;              // theta,phi                                // 2x 32bit float
    float wavelength;        // photon wavelength                  //    32bit float
    float cherenkovDist;     // Cherenkov distance travelled    //    32bit float
    uint32_t numScatters;    // number of scatters                 //    32bit unsigned
    float weight;            //    32bit float
    uint32_t identifier;     //    32bit unsigned
    short stringID;          //    16bit signed
    ushort omID;             //    16bit unsigned
    float4 startPosAndTime;  // 4x 32bit float
    float2 startDir;         // 2x 32bit float
    float groupVelocity;     //    32bit float
    float distInAbsLens;     //    32bit float
                             // total: 20x 32bit float = 80 bytes

    __host__ __device__ I3CLSimPhotonCuda() : posAndTime{0, 0, 0, 0} {}

    __host__ I3CLSimPhotonCuda(const I3CLSimPhoton &i3photon)
        : posAndTime{i3photon.GetPosX(), i3photon.GetPosY(), i3photon.GetPosZ(), i3photon.GetTime()},
          dir{i3photon.GetDirTheta(), i3photon.GetDirPhi()},
          wavelength(i3photon.GetWavelength()),
          cherenkovDist(i3photon.GetCherenkovDist()),
          numScatters(i3photon.GetNumScatters()),
          weight(i3photon.weight),
          identifier(i3photon.GetID()),
          stringID(i3photon.GetStringID()),
          omID(i3photon.GetOMID()),
          startPosAndTime{i3photon.GetStartPosX(), i3photon.GetStartPosY(), i3photon.GetStartPosZ(),
                          i3photon.GetStartTime()},
          startDir{i3photon.GetStartDirTheta(), i3photon.GetStartDirPhi()},
          groupVelocity(i3photon.GetGroupVelocity()),
          distInAbsLens(i3photon.GetDistInAbsLens())
    {
    }
};

#endif  // DATASTRUCCUDA_CUH
