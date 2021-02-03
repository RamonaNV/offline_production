

#include <sstream>

#include "cuda/I3CLSimCUDADevice.h"

#include <icetray/python/list_indexing_suite.hpp>
#include <icetray/python/copy_suite.hpp>

using namespace boost::python;
namespace bp = boost::python;

bool I3CLSimCUDADevice_equalWrap(const I3CLSimCUDADevice &a, const I3CLSimCUDADevice &b)
{
    return a==b;
}

static std::string 
I3CLSimCUDADevice_repr(const I3CLSimCUDADevice& s)
{
    std::ostringstream oss;

    oss << "I3CLSimCUDADevice('" << s.GetDeviceName() << "')";

    return oss.str();
}

static std::string 
I3CLSimCUDADevice_prettyprint(const I3CLSimCUDADevice& s)
{
    std::ostringstream oss;

    oss << "[ I3CLSimCUDADevice name : "  << s.GetDeviceName() << std::endl
        << std::endl
        << "                           cpu : "  << (s.IsCPU()?"YES":"NO") << std::endl
        << "                           gpu : "  << (s.IsGPU()?"YES":"NO") << std::endl
        << "               maxComputeUnits : "  << s.GetMaxComputeUnits() << std::endl
        << "               maxWorkItemSize : "  << s.GetMaxWorkItemSize() << std::endl
        << "              maxWorkGroupSize : "  << s.GetMaxWorkGroupSize() << std::endl
        << "             maxClockFrequency : "  << s.GetMaxClockFrequencyMhz() << "MHz" << std::endl
        << "                 globalMemSize : "  << static_cast<double>(s.GetGlobalMemSize())/1024./1024. << "MiB" << std::endl
        << "         maxConstantBufferSize : "  << static_cast<double>(s.GetMaxConstantBufferSize())/1024. << "kiB" << std::endl
        << "                  localMemSize : "  << static_cast<double>(s.GetLocalMemSize())/1024. << "kiB" << std::endl
        << "             dedicatedLocalMem : "  << (s.HasDedicatedLocalMem()?"YES":"NO") << std::endl
        << "        errorCorrectionSupport : "  << (s.HasErrorCorrectionSupport()?"YES":"NO") << std::endl
        << "                     available : "  << (s.IsAvailable()?"YES":"NO") << std::endl
        << "                        vendor : "  << s.GetVendor() << std::endl
        << "]" ;

    return oss.str();
}

void register_I3CLSimCUDADevice()
{
    {
        bp::scope I3CLSimCUDADevice_scope = 
        bp::class_<I3CLSimCUDADevice, boost::shared_ptr<I3CLSimCUDADevice> >
        ("I3CLSimCUDADevice", bp::no_init)

        .def("__eq__", &I3CLSimCUDADevice_equalWrap)
        .def("__str__", &I3CLSimCUDADevice_prettyprint)
        .def("__repr__", &I3CLSimCUDADevice_repr)

        .def(copy_suite<I3CLSimCUDADevice>())

        .def("GetAllDevices", &I3CLSimCUDADevice::GetAllDevices)
        .staticmethod("GetAllDevices")

        .add_property("device", &I3CLSimCUDADevice::GetDeviceName)

        .def("IsCPU", &I3CLSimCUDADevice::IsCPU)
        .add_property("cpu", &I3CLSimCUDADevice::IsCPU)
        .def("IsGPU", &I3CLSimCUDADevice::IsGPU)
        .add_property("gpu", &I3CLSimCUDADevice::IsGPU)
        .def("GetMaxComputeUnits", &I3CLSimCUDADevice::GetMaxComputeUnits)
        .add_property("maxComputeUnits", &I3CLSimCUDADevice::GetMaxComputeUnits)
        .def("GetMaxWorkItemSize", &I3CLSimCUDADevice::GetMaxWorkItemSize)
        .add_property("maxWorkItemSize", &I3CLSimCUDADevice::GetMaxWorkItemSize)
        .def("GetMaxWorkGroupSize", &I3CLSimCUDADevice::GetMaxWorkGroupSize)
        .add_property("maxWorkGroupSize", &I3CLSimCUDADevice::GetMaxWorkGroupSize)
        .def("GetMaxClockFrequencyMhz", &I3CLSimCUDADevice::GetMaxClockFrequencyMhz)
        .add_property("maxClockFrequencyMhz", &I3CLSimCUDADevice::GetMaxClockFrequencyMhz)
        .def("GetGlobalMemSize", &I3CLSimCUDADevice::GetGlobalMemSize)
        .add_property("globalMemSize", &I3CLSimCUDADevice::GetGlobalMemSize)
        .def("GetMaxConstantBufferSize", &I3CLSimCUDADevice::GetMaxConstantBufferSize)
        .add_property("maxConstantBufferSize", &I3CLSimCUDADevice::GetMaxConstantBufferSize)
        .def("GetLocalMemSize", &I3CLSimCUDADevice::GetLocalMemSize)
        .add_property("localMemSize", &I3CLSimCUDADevice::GetLocalMemSize)
        .def("HasDedicatedLocalMem", &I3CLSimCUDADevice::HasDedicatedLocalMem)
        .add_property("dedicatedLocalMem", &I3CLSimCUDADevice::HasDedicatedLocalMem)
        .def("HasErrorCorrectionSupport", &I3CLSimCUDADevice::HasErrorCorrectionSupport)
        .add_property("errorCorrectionSupport", &I3CLSimCUDADevice::HasErrorCorrectionSupport)
        .def("IsAvailable", &I3CLSimCUDADevice::IsAvailable)
        .add_property("available", &I3CLSimCUDADevice::IsAvailable)
        // .def("GetVendor", &I3CLSimCUDADevice::GetVendor)
        // .add_property("vendor", &I3CLSimCUDADevice::GetVendor)
        // .def("GetDriverVersion", &I3CLSimCUDADevice::GetDriverVersion)
        // .add_property("driverVersion", &I3CLSimCUDADevice::GetDriverVersion)
        // .def("GetDeviceVersion", &I3CLSimCUDADevice::GetDeviceVersion)
        // .add_property("deviceVersion", &I3CLSimCUDADevice::GetDeviceVersion)
        // .def("GetExtensions", &I3CLSimCUDADevice::GetExtensions)
        // .add_property("extensions", &I3CLSimCUDADevice::GetExtensions)
        //
        //
        // .def("GetUseNativeMath", &I3CLSimCUDADevice::GetUseNativeMath)
        // .def("SetUseNativeMath", &I3CLSimCUDADevice::SetUseNativeMath)
        // .add_property("useNativeMath", &I3CLSimCUDADevice::GetExtensions, &I3CLSimCUDADevice::SetUseNativeMath)
        //
        // .def("GetApproximateNumberOfWorkItems", &I3CLSimCUDADevice::GetApproximateNumberOfWorkItems)
        // .def("SetApproximateNumberOfWorkItems", &I3CLSimCUDADevice::SetApproximateNumberOfWorkItems)
        // .add_property("approximateNumberOfWorkItems", &I3CLSimCUDADevice::GetApproximateNumberOfWorkItems, &I3CLSimCUDADevice::SetApproximateNumberOfWorkItems)
        ;
    }
    

#if 1
    class_<I3CLSimCUDADeviceSeries, I3CLSimCUDADeviceSeriesPtr>("I3CLSimCUDADeviceSeries")
    .def(list_indexing_suite<I3CLSimCUDADeviceSeries>())
    ;
    
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimCUDADevice>, boost::shared_ptr<const I3CLSimCUDADevice> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimCUDADeviceSeries>, boost::shared_ptr<const I3CLSimCUDADeviceSeries> >();

    from_python_sequence<I3CLSimCUDADeviceSeries, variable_capacity_policy>();
#endif
}
