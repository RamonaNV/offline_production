
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

#ifndef EMITTERSTEPCUDA_CUH
#define EMITTERSTEPCUDA_CUH

// unused.

#include <I3CLSimStep.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <math.h>
#include <vector_functions.h>

#include <mathstuff.cuh>

// undone class copy of I3CLSimStep to use on host and device.
struct EmitterStepCUDA  // private I3CLSimStep
{
   public:
    __host__ __device__ EmitterStepCUDA()
        : posAndTime({0, 0, 0, 0}),
          dirAndLengthAndBeta{0, 0, 0, 0},
          numPhotons(0),
          weight(1.0),
          identifier(0),
          sourceType(0)
    {
    }

    __host__ EmitterStepCUDA(const I3CLSimStep &i3step)
        : posAndTime{i3step.GetPosX(), i3step.GetPosY(), i3step.GetPosZ(), i3step.GetTime()},
          dirAndLengthAndBeta{i3step.GetDirTheta(), i3step.GetDirPhi(), i3step.GetLength(), i3step.GetBeta()},
          numPhotons(i3step.numPhotons),
          weight(i3step.weight),
          identifier(i3step.identifier),
          sourceType(i3step.sourceType)
    {
    }

    __host__ __device__ __forceinline__ float GetPosX() const { return ((const float *)&posAndTime)[0]; }
    __host__ __device__ __forceinline__ float GetPosY() const { return ((const float *)&posAndTime)[1]; }
    __host__ __device__ __forceinline__ float GetPosZ() const { return ((const float *)&posAndTime)[2]; }
    __host__ __device__ __forceinline__ float GetTime() const { return ((const float *)&posAndTime)[3]; }
    __host__ __device__ __forceinline__ float GetDirTheta() const { return ((const float *)&dirAndLengthAndBeta)[0]; }
    __host__ __device__ __forceinline__ float GetDirPhi() const { return ((const float *)&dirAndLengthAndBeta)[1]; }
    __host__ __device__ __forceinline__ float GetLength() const { return ((const float *)&dirAndLengthAndBeta)[2]; }
    __host__ __device__ __forceinline__ float GetBeta() const { return ((const float *)&dirAndLengthAndBeta)[3]; }
    __host__ __device__ __forceinline__ uint32_t GetNumPhotons() const { return numPhotons; }
    __host__ __device__ __forceinline__ float GetWeight() const { return weight; }
    __host__ __device__ __forceinline__ uint32_t GetID() const { return identifier; }
    __host__ __device__ __forceinline__ uint8_t GetSourceType() const { return sourceType; }

    __host__ __forceinline__ I3PositionPtr GetPos() const
    {
        return I3PositionPtr(new I3Position(((const float *)&posAndTime)[0], ((const float *)&posAndTime)[1],
                                            ((const float *)&posAndTime)[2]));
    }

    __host__ __forceinline__ I3DirectionPtr GetDir() const
    {
        I3DirectionPtr retval(new I3Direction());
        retval->SetThetaPhi(((const float *)&dirAndLengthAndBeta)[0], ((const float *)&dirAndLengthAndBeta)[1]);
        return retval;
    }

    __host__ __device__ __forceinline__ void SetPosX(const float &val) { ((float *)&posAndTime)[0] = val; }
    __host__ __device__ __forceinline__ void SetPosY(const float &val) { ((float *)&posAndTime)[1] = val; }
    __host__ __device__ __forceinline__ void SetPosZ(const float &val) { ((float *)&posAndTime)[2] = val; }
    __host__ __device__ __forceinline__ void SetTime(const float &val) { ((float *)&posAndTime)[3] = val; }
    __host__ __device__ __forceinline__ void SetDirTheta(const float &val) { ((float *)&dirAndLengthAndBeta)[0] = val; }
    __host__ __device__ __forceinline__ void SetDirPhi(const float &val) { ((float *)&dirAndLengthAndBeta)[1] = val; }
    __host__ __device__ __forceinline__ void SetLength(const float &val) { ((float *)&dirAndLengthAndBeta)[2] = val; }
    __host__ __device__ __forceinline__ void SetBeta(const float &val) { ((float *)&dirAndLengthAndBeta)[3] = val; }
    __host__ __device__ __forceinline__ void SetNumPhotons(const uint32_t &val) { numPhotons = val; }
    __host__ __device__ __forceinline__ void SetWeight(const float &val) { weight = val; }
    __host__ __device__ __forceinline__ void SetID(const uint32_t &val) { identifier = val; }
    __host__ __device__ __forceinline__ void SetSourceType(const uint8_t &val) { sourceType = val; }

    __host__ __device__ __forceinline__ void SetDir(const double &x, const double &y, const double &z)
    {
        ((float *)&dirAndLengthAndBeta)[0] = CalcTheta(x, y, z);
        ((float *)&dirAndLengthAndBeta)[1] = CalcPhi(x, y, z);
    }

   private:
    float4 posAndTime;           // x,y,z,time
    float4 dirAndLengthAndBeta;  // theta,phi,length,beta
    uint32_t numPhotons;
    float weight;
    uint32_t identifier;
    uint8_t sourceType;
    /*
    private: Todo
        friend class icecube::serialization::access;
        template <class Archive> void load(Archive & ar, unsigned version);
        template <class Archive> void save(Archive & ar, unsigned version) const;
        I3_SERIALIZATION_SPLIT_MEMBER();*/
};

#endif  // EMITTERSTEPCUDA_CUH
