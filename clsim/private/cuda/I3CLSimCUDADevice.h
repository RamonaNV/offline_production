
#ifndef I3CLSIMCUDADEVICE_H_INCLUDED
#define I3CLSIMCUDADEVICE_H_INCLUDED

#include <icetray/I3PointerTypedefs.h>

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>

struct cudaDeviceProp;

/**
 * @brief Describes an OpenCL platform, device and
 * a few device-specific parameters
 */
class I3CLSimCUDADevice
{
public:

    static boost::shared_ptr<std::vector<I3CLSimCUDADevice> > GetAllDevices();

    std::string GetDeviceName() const;

    // device settings
    // inline void SetUseNativeMath(bool useit) {useNativeMath_=useit;}
    // inline bool GetUseNativeMath() const {return useNativeMath_;}
    // inline void SetApproximateNumberOfWorkItems(uint32_t newnumber) {approximateNumberOfWorkItems_=newnumber;}
    // inline uint32_t GetApproximateNumberOfWorkItems() const {return approximateNumberOfWorkItems_;}

    // device properties (from OpenCL)
    bool IsCPU() const { return false; };
    bool IsGPU() const { return true; };
    std::size_t GetMaxComputeUnits() const;
    std::size_t GetMaxWorkItemSize() const;
    std::size_t GetMaxWorkGroupSize() const;
    std::size_t GetMaxClockFrequencyMhz() const;
    std::size_t GetGlobalMemSize() const;
    std::size_t GetMaxConstantBufferSize() const;
    std::size_t GetLocalMemSize() const;
    bool HasDedicatedLocalMem() const;
    bool HasErrorCorrectionSupport() const;
    bool IsAvailable() const;
    std::string GetVendor() const;
    int GetDeviceNumber() const { return deviceNumber_; };

private:
    I3CLSimCUDADevice(int number);
    I3CLSimCUDADevice() = delete; // no default construction

    // uint32_t approximateNumberOfWorkItems_;
    int deviceNumber_;
    boost::shared_ptr<cudaDeviceProp> device_;

    friend bool operator==(const I3CLSimCUDADevice &, const I3CLSimCUDADevice &);
};
bool operator==(const I3CLSimCUDADevice &a, const I3CLSimCUDADevice &b);

typedef std::vector<I3CLSimCUDADevice> I3CLSimCUDADeviceSeries;

I3_POINTER_TYPEDEFS(I3CLSimCUDADevice);
I3_POINTER_TYPEDEFS(I3CLSimCUDADeviceSeries);

#endif //I3CLSIMCUDADEVICE_H_INCLUDED
