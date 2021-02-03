#include <cuda/I3CLSimCUDADevice.h>

I3CLSimCUDADevice::I3CLSimCUDADevice(int number) : deviceNumber_(number), device_(new cudaDeviceProp)
{}

I3CLSimCUDADeviceSeriesPtr I3CLSimCUDADevice::GetAllDevices()
{
    I3CLSimCUDADeviceSeriesPtr devices(new I3CLSimCUDADeviceSeries());
    int count = 0;
    cudaError_t err = cudaSuccess;
    if ((err = cudaGetDeviceCount(&count)) != cudaSuccess)
        throw std::runtime_error(cudaGetErrorString(err));
    for (int i=0; i < count; i++) {
        devices->push_back(I3CLSimCUDADevice(i));
        cudaGetDeviceProperties(devices->back().device_.get(), i);
    }

    return devices;
}

// comparison
bool operator==(const I3CLSimCUDADevice &a, const I3CLSimCUDADevice &b)
{
    return memcmp(&(a.device_->uuid), &(b.device_->uuid), sizeof(cudaUUID_t)) == 0;
}

// device information
std::string I3CLSimCUDADevice::GetDeviceName() const {return device_->name;}
std::size_t I3CLSimCUDADevice::GetMaxComputeUnits() const {return device_->multiProcessorCount;}
std::size_t I3CLSimCUDADevice::GetMaxWorkItemSize() const {return device_->maxGridSize[0];}
std::size_t I3CLSimCUDADevice::GetMaxWorkGroupSize() const {return device_->maxThreadsPerBlock;}
std::size_t I3CLSimCUDADevice::GetMaxClockFrequencyMhz() const {return device_->clockRate / 1000;}
std::size_t I3CLSimCUDADevice::GetGlobalMemSize() const {return device_->totalGlobalMem;}
std::size_t I3CLSimCUDADevice::GetMaxConstantBufferSize() const {return device_->totalConstMem;}
std::size_t I3CLSimCUDADevice::GetLocalMemSize() const {return device_->sharedMemPerBlock;}
// maybe?
bool I3CLSimCUDADevice::HasDedicatedLocalMem() const {return true;}
bool I3CLSimCUDADevice::HasErrorCorrectionSupport() const {return device_->ECCEnabled;}
// also maybe?
bool I3CLSimCUDADevice::IsAvailable() const {return true;}
std::string I3CLSimCUDADevice::GetVendor() const {return "NVIDIA";}

