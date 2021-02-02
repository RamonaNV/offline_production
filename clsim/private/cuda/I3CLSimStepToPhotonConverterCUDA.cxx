/**
 * Copyright (c) 2011, 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
 * OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 *
 * $Id$
 *
 * @file I3CLSimStepToPhotonConverterCUDA.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include "I3CLSimStepToPhotonConverterCUDA.h"

#include <inttypes.h>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <cmath>

#define PRINTLC printf("thread 0 - in line %d \n", __LINE__);
#include <chrono>
#include <ctime>
#include <propagationKernelSource.cuh>
#include <random>

// debugging: show GPUtime/photon
#define DUMP_STATISTICS

#include <icetray/I3Units.h>
#include <stdlib.h>

#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <limits>
#include <sstream>
#include <string>

#include "opencl/mwcrng_init.h"

#include "cudaBenchmarkSettings.h"

const bool I3CLSimStepToPhotonConverterCUDA::default_useNativeMath = true;

I3CLSimStepToPhotonConverterCUDA::I3CLSimStepToPhotonConverterCUDA(I3RandomServicePtr randomService,
                                                                       bool useNativeMath)
    : openCLStarted_(false),
      queueToOpenCL_(new I3CLSimQueue<ToOpenCLPair_t>(2)),
      queueFromOpenCL_(new I3CLSimQueue<I3CLSimStepToPhotonConverter::ConversionResult_t>(2)),
      randomService_(randomService),
      initialized_(false),
      compiled_(false),
      useNativeMath_(useNativeMath),
      disableDoubleBuffering_(true),
      doublePrecision_(false),
      pancakeFactor_(1.),
      maxWorkgroupSize_(0),
      workgroupSize_(0),
      maxNumWorkitems_(10240)
{
    if (!randomService_) log_fatal("You need to supply a I3RandomService.");
}

I3CLSimStepToPhotonConverterCUDA::~I3CLSimStepToPhotonConverterCUDA()
{
    if (openCLThreadObj_) {
        if (openCLThreadObj_->joinable()) {
            log_debug("Stopping the OpenCL worker thread..");

            openCLThreadObj_->interrupt();

            openCLThreadObj_->join();  // wait for it indefinitely

            log_debug("OpenCL worker thread stopped.");
        }

        openCLThreadObj_.reset();
    }
}

uint64_t I3CLSimStepToPhotonConverterCUDA::GetMaxWorkgroupSize() const
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterCUDA already initialized!");

    if (!compiled_)
        throw I3CLSimStepToPhotonConverter_exception("You need to compile the kernel first. Call Compile().");

    return maxWorkgroupSize_;
}

void I3CLSimStepToPhotonConverterCUDA::SetDevice(const I3CLSimCUDADevice &device)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterCUDA already initialized!");

    device_ = I3CLSimCUDADevicePtr(new I3CLSimCUDADevice(device));  // make a copy of the device
}

void I3CLSimStepToPhotonConverterCUDA::SetWorkgroupSize(std::size_t val)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterCUDA already initialized!");

    workgroupSize_ = val;
}

void I3CLSimStepToPhotonConverterCUDA::SetMaxNumWorkitems(std::size_t val)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterCUDA already initialized!");

    if (val <= 0) throw I3CLSimStepToPhotonConverter_exception("Invalid maximum number of work items!");

    maxNumWorkitems_ = val;
}

std::size_t I3CLSimStepToPhotonConverterCUDA::GetWorkgroupSize() const
{
    if (workgroupSize_ == 0) {
        if (!compiled_)
            throw I3CLSimStepToPhotonConverter_exception(
                "Automatic workgroup size cannot be returned before Compile() has been called!");

        return maxWorkgroupSize_;
    }

    return workgroupSize_;
}

std::size_t I3CLSimStepToPhotonConverterCUDA::GetMaxNumWorkitems() const { return maxNumWorkitems_; }

void I3CLSimStepToPhotonConverterCUDA::Initialize()
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterCUDA already initialized!");

    // use the maximum workgroup size (==granularity) if none has been configured
    // FIXME value picked out of a hat; what's the equivalent of getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE> ?
    maxWorkgroupSize_ = 512*4;
    if (workgroupSize_ == 0) workgroupSize_ = maxWorkgroupSize_;

    if (workgroupSize_ > maxWorkgroupSize_) throw I3CLSimStepToPhotonConverter_exception("Workgroup size too large!");

    const unsigned int numBuffers = disableDoubleBuffering_ ? 1 : 2;

    log_warn("setting buffer sizes");
    // For GPUs, choose the largest number of work items such that they still
    // fit in device memory
    if (device_->IsGPU()) {
        size_t sizePerWorkitem =
            numBuffers * (sizeof(I3CLSimStep)     // the input step
                          + 2 * sizeof(uint64_t)  // MWC multipliers
                          + 10u*sizeof(I3CLSimPhoton)  // the output buffer
                          );
        maxNumWorkitems_ =
            (device_->GetGlobalMemSize() - numBuffers * sizeof(uint32_t)) / sizePerWorkitem;
        size_t numMultipliers = 6139850;
        if (maxNumWorkitems_ > numMultipliers) {
            log_info_stream("Limiting number of work items to " << numMultipliers
                                                                << " (maximum number of prime multipliers)");
            maxNumWorkitems_ = numMultipliers;
        }
        // Choose a bunch size that is a multiple of both the number of cores
        // and the workgroup size
        size_t granularity = device_->GetMaxComputeUnits() * workgroupSize_;
        maxNumWorkitems_ = (maxNumWorkitems_ / granularity) * granularity;
    }

    if (maxNumWorkitems_ % workgroupSize_ != 0)
        throw I3CLSimStepToPhotonConverter_exception(
            "The maximum number of work items (" + boost::lexical_cast<std::string>(maxNumWorkitems_) +
            ") must be a multiple of the workgroup size (" + boost::lexical_cast<std::string>(workgroupSize_) + ").");

    // start with a maximum number of output photons of the same size as the number of
    // input steps. Should be plenty..
    maxNumOutputPhotons_ = static_cast<uint32_t>(std::min(maxNumWorkitems_*10, static_cast<std::size_t>(std::numeric_limits<uint32_t>::max())));

    log_debug("basic OpenCL setup done.");

    // set up rng
    log_warn("Setting up RNG for %zu workitems.", maxNumWorkitems_);

    MWC_RNG_x.resize(maxNumWorkitems_);
    MWC_RNG_a.resize(maxNumWorkitems_);

    if (init_MWC_RNG(&(MWC_RNG_x[0]), &(MWC_RNG_a[0]), maxNumWorkitems_, randomService_) != 0)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterCUDA already initialized!");

    log_debug("RNG is set up..");

    log_debug("Starting the OpenCL worker thread..");
    openCLStarted_ = false;

    openCLThreadObj_ = boost::shared_ptr<boost::thread>(
        new boost::thread(boost::bind(&I3CLSimStepToPhotonConverterCUDA::OpenCLThread, this)));

    // wait for startup
    {
        boost::unique_lock<boost::mutex> guard(openCLStarted_mutex_);
        for (;;) {
            if (openCLStarted_) break;
            openCLStarted_cond_.wait(guard);
        }
    }

    log_debug("CUDA worker thread started.");

    log_debug("CUDA setup complete.");

    initialized_ = true;
}

void I3CLSimStepToPhotonConverterCUDA::OpenCLThread()
{
    // do not interrupt this thread by default
    boost::this_thread::disable_interruption di;

    try {
        //   OpenCLThread_impl(di);
        ThreadyThread(di);
        PRINTLC
    } catch (...) {  // any exceptions?
        std::cerr << "OpenCL worker thread died unexpectedly.." << std::endl;
        exit(0);  // get out as quickly as possible, we probably just had a FATAL error anyway..
        throw;    // will never be reached
    }
}

bool I3CLSimStepToPhotonConverterCUDA::IsInitialized() const { return initialized_; }

void I3CLSimStepToPhotonConverterCUDA::SetEnableDoubleBuffering(bool value)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterCUDA already initialized!");

    disableDoubleBuffering_ = (!value);
}

bool I3CLSimStepToPhotonConverterCUDA::GetEnableDoubleBuffering() const { return (!disableDoubleBuffering_); }

void I3CLSimStepToPhotonConverterCUDA::SetDoublePrecision(bool value)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterCUDA already initialized!");

    doublePrecision_ = value;
}

bool I3CLSimStepToPhotonConverterCUDA::GetDoublePrecision() const { return doublePrecision_; }

void I3CLSimStepToPhotonConverterCUDA::SetDOMPancakeFactor(double value)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterCUDA already initialized!");

    pancakeFactor_ = value;
}

double I3CLSimStepToPhotonConverterCUDA::GetDOMPancakeFactor() const { return pancakeFactor_; }

void I3CLSimStepToPhotonConverterCUDA::SetWlenGenerators(
    const std::vector<I3CLSimRandomValueConstPtr> &wlenGenerators)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterCUDA already initialized!");

    wlenGenerators_ = wlenGenerators;
}

void I3CLSimStepToPhotonConverterCUDA::SetWlenBias(I3CLSimFunctionConstPtr wlenBias)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterCUDA already initialized!");

    wlenBias_ = wlenBias;
}

void I3CLSimStepToPhotonConverterCUDA::SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterCUDA already initialized!");

    mediumProperties_ = mediumProperties;
}

void I3CLSimStepToPhotonConverterCUDA::SetGeometry(I3CLSimSimpleGeometryConstPtr geometry)
{
    if (initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterCUDA already initialized!");
    if (!geometry_) {
        return;
    }
    // The kernel has IC86 geometry hard-coded with contiguous, positive
    // string/dom numbering. Fail loudly if this does not match the supplied
    // geometry.
    if (geometry_->size() != 86u*60u) {
        log_fatal_stream("IC86 geometry hard-coded, but you configured " << geometry_->size() << " modules");
    }
    for (unsigned i=0; i < geometry_->size(); i++) {
        if (geometry_->GetStringID(i) != int(i)/60) {
            log_fatal_stream("String " << geometry_->GetStringID(i) << " where " << i/60 << " should be. Is this really IC86?");
        }
        if (geometry_->GetDomID(i) != i % 60) {
            log_fatal_stream("OM " << geometry_->GetDomID(i) << " where " << i%60 << " should be. Is this really IC86?");
        }
    }

    geometry_ = geometry;
}



void I3CLSimStepToPhotonConverterCUDA::EnqueueSteps(I3CLSimStepSeriesConstPtr steps, uint32_t identifier)
{
    if (!initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterCUDA is not initialized!");

    if (!steps) throw I3CLSimStepToPhotonConverter_exception("Steps pointer is (null)!");

    if (steps->empty()) throw I3CLSimStepToPhotonConverter_exception("Steps are empty!");

    if (steps->size() > maxNumWorkitems_)
        throw I3CLSimStepToPhotonConverter_exception("Number of steps is greater than maximum number of work items!");

    if (steps->size() % workgroupSize_ != 0)
        throw I3CLSimStepToPhotonConverter_exception("The number of steps is not a multiple of the workgroup size!");

    queueToOpenCL_->Put(make_pair(identifier, steps));
}

std::size_t I3CLSimStepToPhotonConverterCUDA::QueueSize() const
{
    if (!initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterCUDA is not initialized!");

    return queueToOpenCL_->size();
}

bool I3CLSimStepToPhotonConverterCUDA::MorePhotonsAvailable() const
{
    if (!initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterCUDA is not initialized!");

    return (!queueFromOpenCL_->empty());
}

// helper
namespace
{
inline void ReplaceStringDOMIndexWithStringDOMIDs(I3CLSimPhotonSeries &photons)
{
    BOOST_FOREACH (I3CLSimPhoton &photon, photons) {
        const int16_t stringIndex = photon.stringID;
        const uint16_t DOMIndex = photon.omID;

        // naive mapping, since IC86 geometry is hard-coded in the CUDA kernel
        const int stringID = stringIndex+1;
        const unsigned int domID = DOMIndex+1;

        if ((stringID < std::numeric_limits<int16_t>::min()) || (stringID > std::numeric_limits<int16_t>::max()))
            log_fatal(
                "Your detector I3Geometry uses a string ID \"%i\". Large IDs like that are currently not supported by "
                "clsim.",
                stringID);

        if (domID > std::numeric_limits<uint16_t>::max())
            log_fatal(
                "Your detector I3Geometry uses a OM ID \"%u\". Large IDs like that are currently not supported by "
                "clsim.",
                domID);

        photon.stringID = static_cast<int16_t>(stringID);
        photon.omID = static_cast<uint16_t>(domID);

        log_trace("Replaced ID (%" PRIi16 "/%" PRIu16 ") with ID (%" PRIi16 "/%" PRIu16 ") (photon @ pos=(%g,%g,%g))",
                  stringIndex, DOMIndex, photon.stringID, photon.omID, photon.GetPosX(), photon.GetPosY(),
                  photon.GetPosZ());
    }
}

}  // namespace

I3CLSimStepToPhotonConverter::ConversionResult_t I3CLSimStepToPhotonConverterCUDA::GetConversionResult()
{
    if (!initialized_)
        throw I3CLSimStepToPhotonConverter_exception("I3CLSimStepToPhotonConverterCUDA is not initialized!");

    ConversionResult_t result = queueFromOpenCL_->Get();

    if (result.photons) {
        ReplaceStringDOMIndexWithStringDOMIDs(*result.photons);
    }

    return result;
}

std::map<std::string, double> I3CLSimStepToPhotonConverterCUDA::GetStatistics() const
{
    std::map<std::string, double> summary;
    boost::unique_lock<boost::mutex> guard(statistics_.mutex);

    const double totalNumPhotonsGenerated = statistics_.total_num_photons_generated;
    const double totalDeviceTime = static_cast<double>(statistics_.total_device_duration) * I3Units::ns;
    const double totalHostTime = static_cast<double>(statistics_.total_host_duration) * I3Units::ns;

    summary["TotalDeviceTime"] = totalDeviceTime;
    summary["TotalHostTime"] = totalHostTime;
    summary["TotalQueueTime"] = statistics_.total_queue_duration;

    summary["DeviceTimePerKernelMean"] = statistics_.device_duration.mean();
    summary["DeviceTimePerKernelStd"] = std::sqrt(statistics_.device_duration.variance());
    summary["HostTimePerKernelMean"] = statistics_.host_duration.mean();
    summary["HostTimePerKernelStd"] = std::sqrt(statistics_.host_duration.variance());
    summary["InputTimePerKernelMean"] = statistics_.input_wait.mean();
    summary["InputTimePerKernelStd"] = std::sqrt(statistics_.input_wait.variance());
    summary["OutputTimePerKernelMean"] = statistics_.output_wait.mean();
    summary["OutputTimePerKernelStd"] = std::sqrt(statistics_.output_wait.variance());

    summary["NumKernelCalls"] = statistics_.total_kernel_calls;
    summary["TotalNumPhotonsGenerated"] = totalNumPhotonsGenerated;
    summary["TotalNumPhotonsAtDOMs"] = statistics_.total_num_photons_atDOMs;

    summary["AverageDeviceTimePerPhoton"] = totalDeviceTime / totalNumPhotonsGenerated;
    summary["AverageHostTimePerPhoton"] = totalHostTime / totalNumPhotonsGenerated;
    summary["DeviceUtilization"] = totalDeviceTime / totalHostTime;

    return summary;
}

void I3CLSimStepToPhotonConverterCUDA::ThreadyThread(boost::this_thread::disable_interruption &di)
{
    uint32_t stepsIdentifier = 0;
    I3CLSimStepSeriesConstPtr steps;

    Kernel kernel(device_->GetDeviceNumber(), maxNumWorkitems_, maxNumOutputPhotons_, MWC_RNG_x, MWC_RNG_a);

    // notify the main thread that everything is set up
    {
        boost::unique_lock<boost::mutex> guard(openCLStarted_mutex_);
        openCLStarted_=true;
    }
    openCLStarted_cond_.notify_all();

    while (true) {
        try {
            boost::this_thread::restore_interruption ri(di);
            log_debug_stream("Waiting for steps");

            std::tie(stepsIdentifier, steps) = queueToOpenCL_->Get();
        } catch (boost::thread_interrupted &i) {
            break;
        }

        i3_assert(steps);
        kernel.uploadSteps(*steps);
        log_debug_stream("Uploaded " << steps->size() << " steps");
        kernel.execute();
        auto photons = kernel.downloadPhotons();
        log_debug_stream("Got " << photons.size() << " photons");

        try {
            boost::this_thread::restore_interruption ri(di);
            queueFromOpenCL_->Put(ConversionResult_t(stepsIdentifier, boost::make_shared<I3CLSimPhotonSeries>(std::move(photons)), nullptr));
        } catch (boost::thread_interrupted &i) {
            break;
        }
    }
}

