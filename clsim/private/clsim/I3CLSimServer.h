
#ifndef CLSIM_I3CLSIMSERVER_H_INCLUDED
#define CLSIM_I3CLSIMSERVER_H_INCLUDED

#include "clsim/I3CLSimStep.h"
#include "clsim/I3CLSimStepToPhotonConverter.h"

class I3CLSimServer {
public:
    I3CLSimServer(const std::string &address, const std::vector<I3CLSimStepToPhotonConverterPtr> &converters);
    ~I3CLSimServer();
    I3CLSimServer(I3CLSimServer &&);
    I3CLSimServer(const I3CLSimServer &) = delete;
    I3CLSimServer& operator=(I3CLSimServer &&);
    I3CLSimServer& operator=(const I3CLSimServer &) = delete;

    std::map<std::string, double> GetStatistics() const;

    /// Get the address the server is bound to. This may differ from the
    /// address passed to the ctor if that address contains wildcards
    /// instructing the system to assign a random port.
    /// See: http://api.zeromq.org/4-3:zmq-tcp
    std::string GetAddress() const;

private:
    class impl;
    std::unique_ptr<impl> impl_;
};

/// @brief Client for communicating with I3CLSimServer
///
/// None of the methods are thread-safe, but EnqueueSteps/EnqueueBarrier
/// and GetConversionResultWithBarrierInfo may be called from different threads
/// to feed and drain the client asynchronously.
class I3CLSimClient {
public:
    I3CLSimClient(const std::string &server_address);
    ~I3CLSimClient();
    I3CLSimClient(I3CLSimClient &&);
    I3CLSimClient(const I3CLSimClient &) = delete;
    I3CLSimClient& operator=(I3CLSimClient &&);
    I3CLSimClient& operator=(const I3CLSimClient &) = delete;

    /// @brief Submit steps for propagation
    ///
    /// This function will block until the server is ready to handle more steps
    void EnqueueSteps(I3CLSimStepSeriesConstPtr steps, uint32_t identifier);
    /// @brief Flush any steps held by the server
    ///
    /// This function will block until the server is ready to handle more steps
    void EnqueueBarrier();
    /// @brief Retrieve propagated photons from the next step bunch
    ///
    /// Bunches will be returned in the order they were enqueued.
    /// barrierWasJustReset will be set to `true` when the last bunch enqueued
    /// before the barrier is returned.
    ///
    /// This function will block until results are available
    I3CLSimStepToPhotonConverter::ConversionResult_t GetConversionResultWithBarrierInfo(bool &barrierWasJustReset);

    std::size_t GetWorkgroupSize() const;
    std::size_t GetMaxNumWorkitems() const;

private:
    class impl;
    std::unique_ptr<impl> impl_;
};

#endif // CLSIM_I3CLSIMSERVER_H_INCLUDED

