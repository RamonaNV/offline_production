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
 * $Id: I3CLSimStepToPhotonConverterCUDA.h 177703 2019-12-05 20:04:47Z jvansanten $
 *
 * @file I3CLSimStepToPhotonConverterCUDA.h
 * @version $Revision: 177703 $
 * @date $Date: 2019-12-05 13:04:47 -0700 (Thu, 05 Dec 2019) $
 * @author Claudio Kopper
 */

#ifndef I3CLSIMSTEPTOPHOTONCONVERTERCUDA_H_INCLUDED
#define I3CLSIMSTEPTOPHOTONCONVERTERCUDA_H_INCLUDED

#include "clsim/I3CLSimStepToPhotonConverter.h"

#include "phys-services/I3RandomService.h"

#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>

#include "clsim/I3CLSimQueue.h"

#include "I3CLSimCUDADevice.h"

#include <vector>
#include <map>
#include <string>
#include <stdexcept>

I3_FORWARD_DECLARATION(I3CLSimFunctionFromTable);
I3_FORWARD_DECLARATION(I3CLSimRandomValueInterpolatedDistribution);

/**
 * @brief Creates photons from a given list of steps and propagates
 * them to a DOM using an OpenCL-enabled algorithm
 */
struct I3CLSimStepToPhotonConverterCUDA : public I3CLSimStepToPhotonConverter
{
public:
    static const bool default_useNativeMath;
    
    I3CLSimStepToPhotonConverterCUDA(I3RandomServicePtr randomService,
                                       bool useNativeMath=default_useNativeMath);
    virtual ~I3CLSimStepToPhotonConverterCUDA();

    /**
     * Sets the workgroup size. A value of 0 
     * uses the maximum possible workgroup size
     * for the kernel.
     *
     * Will throw if already initialized.
     */
    void SetWorkgroupSize(std::size_t val);

    /**
     * Sets the number of parallel work items.
     *
     * Will throw if already initialized.
     */
    void SetMaxNumWorkitems(std::size_t val);

    /**
     * Gets the current workgroup size.
     */
    std::size_t GetWorkgroupSize() const;
    
    /**
     * Gets the number of parallel work items.
     */
    std::size_t GetMaxNumWorkitems() const;
    
    /**
     * Sets the CUDA device.
     *
     * Will throw if already initialized.
     */
    void SetDevice(const I3CLSimCUDADevice &device);
    
    /**
     * Returns the maximum workgroup size for the
     * current kernel.
     *
     * Will throw if compilation fails.
     * Will throw if initialized.
     */
    uint64_t GetMaxWorkgroupSize() const;

    /**
     * Disables or enables double-buffered
     * GPU usage. Double buffering will use
     * two command queues and two sets of input
     * and output buffers in order to transfer
     * data to the GPU while a kernel is executing
     * on the other buffer.
     *
     * This has been observed to yield empty results
     * results on older drivers for the nVidia
     * architecture, so it is disabled by default.
     *
     * Before enabling this for a certain driver/hardware
     * combination, make sure that both correct results
     * are returned. Most of the time the second buffer
     * results are always empty, so this error should be
     * easy to observe.
     *
     * Will throw if already initialized.
     */
    void SetEnableDoubleBuffering(bool value);

    /**
     * Returns true if double buffering is enabled.
     */
    bool GetEnableDoubleBuffering() const;

    /**
     * Enables double-precision support in the
     * kernel. This slows down calculations and
     * requires more memory.
     *
     * The performance hit is minimal on CPUs
     * but up to an order of magnitude on GPUs.
     *
     * Will throw if already initialized.
     */
    void SetDoublePrecision(bool value);
        
    /**
     * Returns true if double precision is enabled.
     */
    bool GetDoublePrecision() const;

    /**
     * Sets the "pancake" factor for DOMs. For standard
     * oversized-DOM simulations, this should be the
     * radius oversizing factor. This will flatten the
     * DOM in the direction parallel to the photon.
     * The DOM will have a pancake-like shape, elongated
     * in the directions perpendicular to the photon direction.
     *
     * The DOM radius (supplied by the geometry) must also include
     * the oversizing factor.
     *
     * Will throw if already initialized.
     */
    void SetDOMPancakeFactor(double value);
    
    /**
     * Returns the "pancake" factor for DOMs.
     */
    double GetDOMPancakeFactor() const;

    /**
     * Sets the wavelength generators. 
     * The first generator (index 0) is assumed to return a Cherenkov
     * spectrum that may have a bias applied to it. This bias factor
     * needs to be set using SetWlenBias().
     * All other generator indices are assumed to be for flasher/laser
     * light generation. During generation, no Cherenkov angle
     * rotation will be applied to those photons with indices >= 1.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetWlenGenerators(const std::vector<I3CLSimRandomValueConstPtr> &wlenGenerators);
    
    /**
     * Sets the wavelength weights. Set this to a constant value
     * of 1 if you do not need biased photon generation.
     * The wavelength spectrum set with SetWlenGenerator()
     * is assumed to have a biassing factor already applied to it.
     * This call sets this factor in order to be able to assign
     * correct weights.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetWlenBias(I3CLSimFunctionConstPtr wlenBias);

    /**
     * Sets the medium properties.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties);
    
    /**
     * Sets the geometry.
     * Will throw if used after the call to Initialize().
     */
    virtual void SetGeometry(I3CLSimSimpleGeometryConstPtr geometry);

    /**
     * Initializes the simulation.
     * Will throw if already initialized.
     */
    virtual void Initialize();

    /**
     * Returns true if initialized.
     * Never throws.
     */
    virtual bool IsInitialized() const;
    
    /**
     * Adds a new I3CLSimStepSeries to the queue.
     * The resulting I3CLSimPhotonSeries can be retrieved from the
     * I3CLSimStepToPhotonConverter after some processing time.
     *
     * Enqueuing a vector after calling EnqueueBarrier 
     * will throw if not all photons have been retrieved.
     *
     * Will throw if not initialized.
     */
    virtual void EnqueueSteps(I3CLSimStepSeriesConstPtr steps, uint32_t identifier);

    /**
     * Reports the current queue size. The queue works asynchronously,
     * so this value will probably have changed once you use it.
     *
     * Will throw if not initialized.
     */
    virtual std::size_t QueueSize() const; 

    /**
     * Returns true if more photons are available.
     * If the return value is false, the current simulation is finished
     * and a new step vector may be set.
     * 
     * Will throw if not initialized.
     */
    virtual bool MorePhotonsAvailable() const;

    /**
     * Returns a bunch of photons stored in a vector<I3CLSimPhoton>.
     *
     * The return value is a pair<uint, vector<I3CLSimPhoton> >.
     * The integer is the same identifier as specified in the call
     * to EnqueueSteps().
     *
     * Might block if no photons are available.
     * 
     * Will throw if not initialized.
     */
    virtual I3CLSimStepToPhotonConverter::ConversionResult_t GetConversionResult();
    
    virtual std::map<std::string, double> GetStatistics() const;
    
    inline double GetTotalDeviceTime() const {boost::unique_lock<boost::mutex> guard(statistics_.mutex); return static_cast<double>(statistics_.total_device_duration);}
    inline double GetTotalHostTime() const {boost::unique_lock<boost::mutex> guard(statistics_.mutex); return static_cast<double>(statistics_.total_host_duration);}
    inline double GetTotalQueueTime() const {boost::unique_lock<boost::mutex> guard(statistics_.mutex); return static_cast<double>(statistics_.total_queue_duration);}
    inline uint64_t GetNumKernelCalls() const {boost::unique_lock<boost::mutex> guard(statistics_.mutex); return statistics_.total_kernel_calls;}
    inline uint64_t GetTotalNumPhotonsGenerated() const {boost::unique_lock<boost::mutex> guard(statistics_.mutex); return statistics_.total_num_photons_generated;}
    inline uint64_t GetTotalNumPhotonsAtDOMs() const {boost::unique_lock<boost::mutex> guard(statistics_.mutex); return statistics_.total_num_photons_atDOMs;}

private:
    typedef std::pair<uint32_t, I3CLSimStepSeriesConstPtr> ToOpenCLPair_t;

    void OpenCLThread();
    void ThreadyThread(boost::this_thread::disable_interruption &);

    // Keep a running mean using Welford's online algorithm
    // See: https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    class metric {
    public:
        metric() : count_(0), mean_(0), mean2_(0) {}
        void update(double value) {
            count_++;
            double delta = value - mean_;
            mean_ += delta / count_;
            double delta2 = value - mean_;
            mean2_ += delta*delta2;
        }
        size_t count() const { return count_; }
        double mean() const { return mean_; }
        double variance() const { return (count_ > 1) ? mean2_/(count_-1) : NAN; }
    private:
        size_t count_;
        double mean_, mean2_;
    };
    struct statistics_bundle {
        statistics_bundle() :
            total_host_duration(0), total_device_duration(0), total_queue_duration(0),
            total_kernel_calls(0), total_num_photons_generated(0), total_num_photons_atDOMs(0) {}
        metric input_wait;
        metric output_wait;
        metric host_duration;
        metric device_duration;
        uint64_t total_host_duration;
        uint64_t total_device_duration;
        uint64_t total_queue_duration;
        uint64_t total_kernel_calls;
        uint64_t total_num_photons_generated;
        uint64_t total_num_photons_atDOMs;
        mutable boost::mutex mutex;
    };
    statistics_bundle statistics_;

    boost::shared_ptr<boost::thread> openCLThreadObj_;
    boost::condition_variable_any openCLStarted_cond_;
    boost::mutex openCLStarted_mutex_;
    bool openCLStarted_;
    
    boost::shared_ptr<I3CLSimQueue<ToOpenCLPair_t> > queueToOpenCL_;
    boost::shared_ptr<I3CLSimQueue<I3CLSimStepToPhotonConverter::ConversionResult_t> > queueFromOpenCL_;

    I3RandomServicePtr randomService_;
    
    bool initialized_;
    bool compiled_;
    I3CLSimRandomValueInterpolatedDistributionConstPtr wlenGenerator_;
    I3CLSimFunctionFromTableConstPtr wlenBias_;
    I3CLSimMediumPropertiesConstPtr mediumProperties_;
    I3CLSimSimpleGeometryConstPtr geometry_;
    
    boost::shared_ptr<I3CLSimCUDADevice> device_;
    bool useNativeMath_;
    bool deviceIsSelected_;
    
    bool disableDoubleBuffering_;
    bool doublePrecision_;
    double pancakeFactor_;

    // maximum workgroup size for current kernel
    uint64_t maxWorkgroupSize_;
    
    // configured workgroup size and maximum number of work items
    std::size_t workgroupSize_;
    std::size_t maxNumWorkitems_;

    // rng state per workitem
    std::vector<uint64_t> MWC_RNG_x;
    std::vector<uint32_t> MWC_RNG_a;

    // Size of output photon storage (maximum amount of photons per step bunch)
    uint32_t maxNumOutputPhotons_;
    
    SET_LOGGER("I3CLSimStepToPhotonConverterCUDA");
};

I3_POINTER_TYPEDEFS(I3CLSimStepToPhotonConverterCUDA);

#endif //I3CLSIMSTEPTOPHOTONCONVERTERCUDA_H_INCLUDED
