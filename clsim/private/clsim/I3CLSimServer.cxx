
#include "icetray/I3Logging.h"
#include "clsim/I3CLSimServer.h"

#ifdef HAS_ZMQ

#include <chrono>
#include <thread>
#include <future>
#include <memory>
#include <map>
#include <list>
#include <queue>
#include <zmq.hpp>

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/integer/common_factor_rt.hpp>
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>

namespace {

std::tuple<zmq::message_t, std::list<zmq::message_t> >
read_message(zmq::socket_t &socket)
{
    zmq::message_t address;
    std::list<zmq::message_t> body;
    
    socket.recv(&address);
    if (!address.more()) {
        body.push_back(std::move(address));
        address.rebuild();
        return std::make_tuple(std::move(address), std::move(body));
    }
    
    do {
        zmq::message_t part;
        socket.recv(&part);
        body.push_back(std::move(part));
    } while ((--body.end())->more());
    
    // remove the delimiter added by REQ or REP sockets
    if (body.front().size() == 0)
        body.pop_front();
        
    return std::make_tuple(std::move(address), std::move(body));
}

// NB: for simple struct-like types such as I3CLSimStep and I3CLSimPhoton, serializing
// to a properly-sized vector<char> is only a few percent slower than a straight memcpy().
// The extra copy to the output message_t 
template<typename T>
zmq::message_t serialize(const T &value)
{
    std::vector<char> output;
    {
        boost::iostreams::filtering_ostream ostream;
        ostream.push(boost::iostreams::back_inserter(output));
        icecube::archive::portable_binary_oarchive poa(ostream);

        poa << value;
    }

    return zmq::message_t(output.data(), output.size());
}

template <typename T>
T deserialize(const zmq::message_t &message)
{
    T value;

    boost::iostreams::filtering_istream istream;
    istream.push(boost::iostreams::array_source(message.data<char>(), message.size()));
    icecube::archive::portable_binary_iarchive pia(istream);

    pia >> value;
    i3_assert( istream.tellg() == static_cast<unsigned>(message.size()) );

    return value;
}

}

namespace I3CLSimServerUtils {

static const zmq::message_t greeting("servus", 6);
static const zmq::message_t barrier("adieu", 5);

};

class I3CLSimServer::impl {
public:
    impl(const std::string &address, const std::vector<I3CLSimStepToPhotonConverterPtr> &converters);
    ~impl();

    std::map<std::string, double> GetStatistics() const;

    std::string GetAddress() const;

private:
    std::vector<I3CLSimStepToPhotonConverterPtr> converters_;
    std::size_t workgroupSize_, maxBunchSize_;
    
    zmq::context_t context_;
    zmq::socket_t control_, heartbeat_;

    SET_LOGGER("I3CLSimServer");

    void ServerThread(const std::string &bindAddress);
    void WorkerThread(unsigned index, unsigned buffer_id);

    std::thread serverThread_;
    std::vector<std::thread> workerThreads_;
    std::string frontendAddress_;
};

I3CLSimServer::impl::impl(const std::string &address, const std::vector<I3CLSimStepToPhotonConverterPtr> &converters) :
    converters_(converters),
    workgroupSize_(0), maxBunchSize_(0),
    context_(1),
    control_(context_, ZMQ_PUB),
    heartbeat_(context_, ZMQ_REP)
{
    control_.bind("inproc://control");
    heartbeat_.bind("inproc://heartbeat");
    
    if (converters_.empty())
        log_fatal("Need at least 1 I3CLSimStepToPhotonConverter");

    // Harmonize bunch sizes
    for (auto &converter : converters_) {
        if (!converter || !converter->IsInitialized())
            log_fatal("All I3CLSimStepToPhotonConverters must be initialized");
        
        if (workgroupSize_ == 0) {
            workgroupSize_ = converter->GetWorkgroupSize();
        } else {
            workgroupSize_ = boost::integer::lcm(workgroupSize_, converter->GetWorkgroupSize());
        }
        
        if (maxBunchSize_ == 0) {
            maxBunchSize_ = converter->GetMaxNumWorkitems();
        } else {
            std::size_t newMaxBunchSize = std::min(maxBunchSize_, converter->GetMaxNumWorkitems());
            std::size_t newMaxBunchSizeWithGranularity = newMaxBunchSize - newMaxBunchSize%workgroupSize_;
            
            if (newMaxBunchSizeWithGranularity == 0)
                log_fatal("maximum bunch sizes are incompatible with kernel work group sizes.");
            maxBunchSize_ = newMaxBunchSizeWithGranularity;
        }
    }
    
    serverThread_ = std::thread(&I3CLSimServer::impl::ServerThread, this, address);
    {
        // Wait for thread to come online so that shutdown messages
        // in the destructor will have an effect
        zmq::message_t ping;
        heartbeat_.recv(&ping);
        frontendAddress_ = deserialize<std::string>(ping);
        heartbeat_.send(ping);
    }
    
    unsigned queueDepth = 1;

    for (unsigned i=0; i < converters_.size(); i++)
        for (unsigned j=0; j < queueDepth; j++) {
            workerThreads_.emplace_back(&I3CLSimServer::impl::WorkerThread, this, i, j);

            zmq::message_t ping;
            heartbeat_.recv(&ping);
            heartbeat_.send(ping);
        }
    
}

I3CLSimServer::impl::~impl()
{
    // send shutdown message to threads
    try {
        control_.send(zmq::message_t());
    } catch (...) {}
    log_debug("Sent shutdown message");
    serverThread_.join();
    for(auto &thread : workerThreads_)
        thread.join();
}

template<typename T>
inline bool send_msg(
  zmq::socket_t &socket,
  const T& payload,
  bool will_continue = false
)
{
  zmq::message_t msg_payload( sizeof(T) );
  std::memcpy(msg_payload.data(), &payload, sizeof(T) );
  return socket.send(msg_payload, will_continue ? ZMQ_SNDMORE : 0); // last message part
}

void I3CLSimServer::impl::ServerThread(const std::string &bindAddress)
{
    zmq::socket_t frontend(context_, ZMQ_ROUTER);
    zmq::socket_t backend(context_, ZMQ_ROUTER);
    {
        // Leave enough room in queues for one bunch of messages
        int bunchBlock = converters_.size()*maxBunchSize_/workgroupSize_;
        frontend.setsockopt(ZMQ_RCVHWM, bunchBlock);
        frontend.setsockopt(ZMQ_SNDHWM, bunchBlock);
        // frontend.setsockopt(ZMQ_SNDHWM, 0);

        frontend.setsockopt(ZMQ_IMMEDIATE, 1);
        frontend.setsockopt(ZMQ_ROUTER_MANDATORY, 1);
        frontend.bind(bindAddress);

        backend.setsockopt(ZMQ_RCVHWM, bunchBlock);
        backend.bind("inproc://worker");
    }

    // Listen for control messages
    zmq::socket_t control(context_, ZMQ_SUB);
    control.connect("inproc://control");
    control.setsockopt(ZMQ_SUBSCRIBE, "", 0);
    {
        // Send bound address to the main thread
        std::string boundAddress;
        size_t size = std::max(size_t(1024), bindAddress.size());
        boundAddress.resize(size);
        frontend.getsockopt(ZMQ_LAST_ENDPOINT, &boundAddress.front(), &size);
        boundAddress.resize(size);

        zmq::socket_t heartbeat(context_, ZMQ_REQ);
        heartbeat.connect("inproc://heartbeat");

        zmq::message_t ping;
        heartbeat.send(serialize(boundAddress));
        heartbeat.recv(&ping);
    }
    
    zmq::pollitem_t items[] = {
      { static_cast<void *>(frontend), 0, ZMQ_POLLIN, 0 },
      { static_cast<void *>(backend),  0, ZMQ_POLLIN, 0 },
      { static_cast<void *>(control),   0, ZMQ_POLLIN, 0 },
    };
    zmq::pollitem_t outitems[] = {
      { static_cast<void *>(frontend), 0, ZMQ_POLLOUT, 0 }
    };
    auto pollout = [](zmq::pollitem_t *item, long timeout=0) {
        zmq::poll(item, 1, timeout);
        return item->revents & ZMQ_POLLOUT;
    };
    
    /// Addresses of idle workers
    std::queue<zmq::message_t > workers;
    /// Map from internal task ID to address of originating client and external ID
    std::map<uint32_t, std::array<zmq::message_t, 2> > clients;

    auto handle_worker = [&]() {
        zmq::message_t address;
        std::list<zmq::message_t> body;
        std::tie(address, body) = read_message(backend);
        
        // A "ready" message; no payload
        if (body.size() == 1 && body.front().size() == 0) {
            workers.push(std::move(address));
            return;
        }
        
        // Otherwise, look up address of originating client
        // and return the payload
        uint32_t client_id = deserialize<uint32_t>(body.back());
        body.pop_back();
        auto destination = clients.find(client_id);
        if (destination == clients.end()) {
            log_error_stream("Unknown client ID " << client_id);
            return;
        }
        log_trace_stream("returning request "<<deserialize<uint32_t>(destination->second[1]));

        try {
            frontend.send(destination->second[0], ZMQ_SNDMORE);
            frontend.send(zmq::message_t(), ZMQ_SNDMORE);
            for (auto &packet : body)
                frontend.send(packet, ZMQ_SNDMORE);
            frontend.send(destination->second[1], 0);
        } catch (zmq::error_t &err) {
            log_error_stream("returning request "<<deserialize<uint32_t>(destination->second[1])<<": "<<err.what());
        }

        clients.erase(destination);
    };

    log_trace_stream("Server thread started on "<<bindAddress);
    while (true) {
        
        // only poll for client messages if a worker is available
        zmq::poll(&items[workers.empty() ? 1 : 0], workers.empty() ? 2 : 3, -1);

        // Message from client
        // Process only if there are available workers
        if (!workers.empty() && (items[0].revents & ZMQ_POLLIN)) {
            
            zmq::message_t address;
            std::list<zmq::message_t> body;
            std::tie(address, body) = read_message(frontend);
            log_trace("got message on frontend");
            
            if ( body.size() == 1 ) {
                if ( body.front() == I3CLSimServerUtils::greeting ) {
                    // a client has just connected. send the desired workgroup batching
                    log_trace_stream("sending batch size");
                    frontend.send(address, ZMQ_SNDMORE);
                    frontend.send(zmq::message_t(), ZMQ_SNDMORE);
                    frontend.send(serialize(std::make_pair(workgroupSize_, maxBunchSize_)), 0);
                } else if ( body.front() == I3CLSimServerUtils::barrier ) {
                    while (workers.size() < workerThreads_.size()) {
                        log_trace_stream(workers.size()<<" workers");
                        handle_worker();
                        log_trace_stream(workers.size()<<" workers");
                    }
                    log_debug_stream("flushing "<<workers.size()<<" workers");

                    while (!workers.empty()) {
                        backend.send(workers.front(), ZMQ_SNDMORE);
                        backend.send(zmq::message_t(), ZMQ_SNDMORE);
                        backend.send(body.front(), 0);
                        workers.pop();
                    }
                    log_debug_stream("flushed");
                } else {
                    log_error_stream("Unknown message "<<std::string(
                        body.front().data<char>(),
                        body.front().data<char>()+body.front().size()));
                }
            } else if ( body.size() == 2 ) {
                // assign an internal ID for later reply to client
                uint32_t client_id(0);
                if (clients.size() > 0)
                    client_id = (--clients.end())->first+1;
                if (clients.find(client_id) != clients.end())
                    log_fatal("Repeated client ID");
                clients.emplace(client_id,
                    std::array<zmq::message_t,2>({{std::move(address), std::move(body.back())}}));
                // drop external ID
                body.pop_back();
            
                // forward to a worker
                log_trace_stream("Forwarding "<<(body.size())<<" packets to worker with id "<<client_id);
                backend.send(workers.front(), ZMQ_SNDMORE);
                backend.send(zmq::message_t(), ZMQ_SNDMORE);
                for (auto &packet : body)
                    backend.send(packet, ZMQ_SNDMORE);
                backend.send(serialize(client_id), 0);
            
                workers.pop();
            }
        }

        // Message from worker
        if ((items[1].revents & ZMQ_POLLIN) && pollout(&outitems[0], 1)) {
            handle_worker();
        }

        // Shutdown signal from main thread
        if (items[2].revents & ZMQ_POLLIN) {
            zmq::message_t dummy;
            control.recv(&dummy);
            log_debug("Server thread shutting down");
            break;
        }
        

    }
}

void I3CLSimServer::impl::WorkerThread(unsigned index, unsigned buffer_id)
{
    // Listen for work
    zmq::socket_t backend(context_, ZMQ_DEALER);
    backend.connect("inproc://worker");
    // Listen for control messages
    zmq::socket_t control(context_, ZMQ_SUB);
    control.connect("inproc://control");
    control.setsockopt(ZMQ_SUBSCRIBE, "", 0);
    // Signal main thread that we're alive
    zmq::socket_t heartbeat(context_, ZMQ_REQ);
    heartbeat.connect("inproc://heartbeat");
    {
        zmq::message_t ping;
        heartbeat.send(ping);
        heartbeat.recv(&ping);
    }
    
    zmq::pollitem_t items[] = {
      { static_cast<void *>(backend), 0, ZMQ_POLLIN, 0 },
      { static_cast<void *>(control), 0, ZMQ_POLLIN, 0 },
    };
    
    // request id <-> bunch id
    boost::bimaps::bimap<
        boost::bimaps::multiset_of<uint32_t>,
        boost::bimaps::multiset_of<uint32_t> > requestsForBunches;
    std::map<uint32_t,
        I3CLSimStepToPhotonConverter::ConversionResult_t> pendingRequests;
    // (request id,step id) <-> local step id
    boost::bimaps::bimap<
        std::pair<uint32_t,uint32_t>, uint32_t> localStepIds;

    // this may wrap, but uniqueness is preserved unless this thread recieves 
    // more than 4 billion unique light sources 
    uint32_t current_bunch_id(0);
    uint32_t current_step_id(0);
    uint32_t pending_bunches(0);
    std::deque<I3CLSimStep> stepStore;

    auto drain = [&](const zmq::message_t &address) {
        I3CLSimStepToPhotonConverter::ConversionResult_t result = 
            converters_[index]->GetConversionResult();
        pending_bunches--;

        // fill results into output
        if (result.photons) {
            auto history = result.photonHistories ?
                result.photonHistories->begin() :
                I3CLSimPhotonHistorySeries::const_iterator();
            for (auto &photon : *result.photons) {
                auto client_id = localStepIds.right.find(photon.GetID());
                // log_trace_stream("photon ID "<<photon.GetID());
                i3_assert( client_id != localStepIds.right.end() );
                auto &request = pendingRequests[client_id->second.first];
                photon.SetID(client_id->second.second);
                request.photons->push_back(photon);
                if (result.photonHistories) {
                    if (!request.photonHistories)
                        request.photonHistories = boost::make_shared<I3CLSimPhotonHistorySeries>();
                    request.photonHistories->push_back(*history++);
                }
            }
        }
        std::deque<uint32_t> finished;
        {
            // Mark this bunch as done
            auto done = requestsForBunches.right.equal_range(result.identifier);
            for (auto bunch=done.first; bunch!=done.second; bunch++) {
                // log_trace_stream("["<<index<<"] request: "<<bunch->second<<" bunch: "<<bunch->first);
                finished.push_back(bunch->second);
            }
            requestsForBunches.right.erase(done.first,done.second);
        }
        for (auto request_id : finished) {
            if (requestsForBunches.left.find(request_id) == requestsForBunches.left.end()) {
                // this request has been fully serviced
                auto request = pendingRequests.find(request_id);
                i3_assert( request != pendingRequests.end() );
                log_trace_stream("["<<index<<"] finished request "<<request_id);
                {
                    zmq::message_t addr;
                    addr.copy(&address);
                    backend.send(addr, ZMQ_SNDMORE);
                    backend.send(zmq::message_t(), ZMQ_SNDMORE);
                }
                if (request->second.photons)
                    backend.send(serialize(request->second.photons), ZMQ_SNDMORE);
                if (request->second.photonHistories)
                    backend.send(serialize(request->second.photonHistories), ZMQ_SNDMORE);
                backend.send(serialize(request->second.identifier), 0);
                // remove cache entry
                pendingRequests.erase(request);
                // remove local id mapping by erasing the range of
                // possible (client_id,step_id) pairs with the given
                // client_id
                localStepIds.left.erase(
                    localStepIds.left.lower_bound(std::make_pair(request_id,std::numeric_limits<uint32_t>::min())),
                    localStepIds.left.upper_bound(std::make_pair(request_id,std::numeric_limits<uint32_t>::max()))
                );
            }
        }
        log_trace_stream("["<<index<<"] finished "<<finished.size()<<" bunches");
    };
    auto submit = [&](I3CLSimStepSeriesPtr bunch, const zmq::message_t &address) {
        // If there are two bunches in the pipeline, wait for one to finish
        if (pending_bunches > 1) {
            drain(address);
        }
        log_trace_stream("["<<index<<"] submitting bunch "<<current_bunch_id);
        converters_[index]->EnqueueSteps(bunch, current_bunch_id++);
        log_trace_stream("["<<index<<"] done! bunch "<<current_bunch_id-1);
        pending_bunches++;
    };

    while (true) {

        // Request work from the server thread
        // (empty message to emulate REQ envelope)
        backend.send(zmq::message_t(),ZMQ_SNDMORE);
        backend.send(zmq::message_t());

        if (int rc = zmq::poll(&items[0], 2, -1) < 0) {
            log_error_stream("ZMQ polling error " << rc);
            break;
        }
        
        // Message from client
        if (items[0].revents & ZMQ_POLLIN) {
            zmq::message_t address;
            std::list<zmq::message_t> body;
            bool flush = false;
            {
                std::tie(address, body) = read_message(backend);
                I3CLSimStepSeriesPtr steps;
                uint32_t request_id(UINT_MAX);
                if ( body.size() == 2 ) {
                    auto msg = body.begin();
                    steps = deserialize<I3CLSimStepSeriesPtr>(*msg++);
                    request_id = deserialize<uint32_t>(*msg++);
                } else {
                    // i3_assert(body.size() == 1);
                }
                
                // Submit job
                if (steps) {
                    for (I3CLSimStep &step: *steps) {
                    
                        auto id = localStepIds.left.find(std::make_pair(request_id,step.GetID()));
                        if (id == localStepIds.left.end()) {
                            id = localStepIds.left.insert(std::make_pair(std::make_pair(
                                request_id,step.GetID()),current_step_id++)).first;
                        }
                        // log_trace_stream("["<<index<<"] request "<<request_id<<" step ID "<<step.GetID()<<" -> "<<id->second);
                        step.SetID(id->second);
                        stepStore.push_back(step);
                    }
                    pendingRequests.emplace(request_id,
                        I3CLSimStepToPhotonConverter::ConversionResult_t(
                            request_id,
                            boost::make_shared<I3CLSimPhotonSeries>()
                        )
                    );
                    // This request will be handled one or more bunches, now
                    while (stepStore.size() >= maxBunchSize_) {
                        // and possibly following ones
                        auto range = std::make_pair(stepStore.begin(),
                            std::next(stepStore.begin(),maxBunchSize_));
                        auto bunch = boost::make_shared<I3CLSimStepSeries>(range.first,range.second);
                        stepStore.erase(range.first,range.second);
                        // log_trace_stream("["<<index<<"] request "<<request_id<<" to bunch "<<current_bunch_id);
                        requestsForBunches.left.insert(std::make_pair(request_id,current_bunch_id));
                        submit(bunch, address);
                    }
                    // This request will be handled in the next bunch
                    if (!stepStore.empty()) {
                        // log_trace_stream("["<<index<<"] request "<<request_id<<" to bunch "<<current_bunch_id);
                        requestsForBunches.left.insert(std::make_pair(request_id,current_bunch_id));
                    }
                } else if (!stepStore.empty()) {
                    log_debug_stream("["<<index<<"] flushing "<<stepStore.size()<<" steps from store");
                    // barrier reached; flush entire step store
                    auto bunch = boost::make_shared<I3CLSimStepSeries>(stepStore.begin(),stepStore.end());
                    stepStore.clear();
                    I3CLSimStep NoOpStepTemplate;
                    NoOpStepTemplate.SetPos(I3Position(0.,0.,0.));
                    NoOpStepTemplate.SetDir(I3Direction(0.,0.,-1.));
                    NoOpStepTemplate.SetTime(0.);
                    NoOpStepTemplate.SetLength(0.);
                    NoOpStepTemplate.SetNumPhotons(0);
                    NoOpStepTemplate.SetWeight(0.);
                    NoOpStepTemplate.SetBeta(1.);
                    size_t padding = (workgroupSize_-(bunch->size()%workgroupSize_))%workgroupSize_;
                    for (size_t i=0; i < padding; i++)
                        bunch->push_back(NoOpStepTemplate);
                    submit(bunch, address);
                    flush = true;
                }
            }
            if (pending_bunches > 0)
                log_trace_stream("["<<index<<"] "<<pending_bunches<<" pending bunches, "<<stepStore.size()<<" steps in store");

            // Get next result (not necessarily from the batches we just enqueued)
            while (pending_bunches > 0 && (flush || converters_[index]->MorePhotonsAvailable())) {
                drain(address);
            }
        }
        
        // Shutdown signal from main thread
        if (items[1].revents & ZMQ_POLLIN) {
            zmq::message_t dummy;
            control.recv(&dummy);
            break;
        }
        
    }
}

std::map<std::string, double> I3CLSimServer::impl::GetStatistics() const
{
    std::map<std::string, double> summary;
    
    for (std::size_t i=0; i<converters_.size(); ++i)
    {
        const std::string postfix = (converters_.size()==1)?"":"_"+boost::lexical_cast<std::string>(i);
        for (auto &v : converters_[i]->GetStatistics()) {
            summary[v.first + postfix] = v.second;
        }
    }
    
    return summary;
}

std::string I3CLSimServer::impl::GetAddress() const
{
    return frontendAddress_;
}

namespace {

SET_LOGGER("I3CLSimClient");

// Since ZMQ sockets are not thread-safe, we use a connector thread to shuffle
// messages between an inbox and and outbox that can each be serviced from
// different threads.
void ClientWorker(zmq::context_t &context, const std::string &serverAddress, int queueDepth=1) {
    zmq::socket_t server(context, ZMQ_DEALER), control(context, ZMQ_SUB),
        inbox(context, ZMQ_PULL), outbox(context, ZMQ_PUSH);
    control.setsockopt(ZMQ_SUBSCRIBE, "", 0);
    inbox.setsockopt(ZMQ_RCVHWM, queueDepth);
    outbox.setsockopt(ZMQ_SNDHWM, queueDepth);
    control.connect("inproc://control");
    inbox.connect("inproc://inbox");
    outbox.connect("inproc://outbox");

    // connect and forward join message to main thread
    server.setsockopt(ZMQ_SNDHWM, queueDepth);
    server.setsockopt(ZMQ_RCVHWM, queueDepth);
    server.connect(serverAddress);
    {
        zmq::socket_t heartbeat(context, ZMQ_REQ);
        heartbeat.connect("inproc://heartbeat");

        zmq::message_t address;
        std::list<zmq::message_t> body;
        std::string greeting = "servus";
        log_trace_stream("connecting to "<<serverAddress);
        server.send(zmq::message_t(greeting.data(), greeting.size()), 0);
        std::tie(address, body) = read_message(server);
        for (auto &msg : body) {
            heartbeat.send(std::move(msg));
        }
        heartbeat.recv(&address);
    }

    zmq::pollitem_t items[] = {
      { static_cast<void *>(server), 0, ZMQ_POLLIN, 0 },
      { static_cast<void *>(inbox),  0, ZMQ_POLLIN, 0 },
      { static_cast<void *>(control),  0, ZMQ_POLLIN, 0 }
    };
    // poll separately for output to avoid busy-waiting when only the output
    // queue is empty
    zmq::pollitem_t outitems[] = {
      { static_cast<void *>(server), 0, ZMQ_POLLOUT, 0 },
      { static_cast<void *>(outbox), 0, ZMQ_POLLOUT, 0 }
    };
    auto pollout = [](zmq::pollitem_t *item, long timeout=0) {
        zmq::poll(item, 1, timeout);
        return item->revents & ZMQ_POLLOUT;
    };
    zmq::message_t lastIdentifier;
    std::set<uint32_t> pendingIdentifiers;
    bool barrierActive = false;
    const zmq::message_t barrier("adieu", 5);
    log_trace_stream("Connector thread running");

    while (true) {

        zmq::poll(&items[0], 3, -1);

        if ((items[0].revents & ZMQ_POLLIN) && pollout(&outitems[1], 1)) {
            zmq::message_t address;
            std::list<zmq::message_t> body;
            std::tie(address, body) = read_message(server);
            i3_assert( body.size() > 0 );
            // check whether the barrier has been reached
            auto id = pendingIdentifiers.find(deserialize<uint32_t>(body.back()));
            if (id == pendingIdentifiers.end()) {
                log_fatal_stream("Unknown identifier");
            } else {
                pendingIdentifiers.erase(id);
            }
            outbox.send(serialize(barrierActive && pendingIdentifiers.empty()), ZMQ_SNDMORE);
            if (pendingIdentifiers.empty()) {
                barrierActive = false;
            }
            for (auto &msg : body) {
                int flags = msg.more() ? ZMQ_SNDMORE : 0;
                outbox.send(std::move(msg), flags);
            }
        }
        
        // Forward a messages to the server only if they will not block
        if ((items[1].revents & ZMQ_POLLIN) && pollout(&outitems[0], 1)) {
            int count(0);
            int flags(0);
            do {
                zmq::message_t msg;
                inbox.recv(&msg);
                if (count == 0 && !msg.more() && msg == I3CLSimServerUtils::barrier) {
                    // barrier encountered
                    barrierActive = true;
                } else if (!msg.more()) {
                    // note last id seen
                    pendingIdentifiers.insert(deserialize<uint32_t>(msg));
                }
                flags = msg.more() ? ZMQ_SNDMORE : 0;
                server.send(std::move(msg), flags);
                count++;
            } while (flags);

            // server will never reply to a flush with no pending steps
            if (barrierActive && pendingIdentifiers.empty()) {
                outbox.send(serialize(barrierActive), ZMQ_SNDMORE);
                outbox.send(serialize(UINT_MAX), 0);
                barrierActive = false;
            }
        }

        if (items[2].revents & ZMQ_POLLIN) {
            zmq::message_t msg;
            control.recv(&msg);
            break;
        }
    }
}

}

class I3CLSimClient::impl {
public:
    impl(const std::string &server_address);
    ~impl();

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

    std::size_t GetWorkgroupSize() const { return workgroupSize_; }
    std::size_t GetMaxNumWorkitems() const { return maxBunchSize_; }
    
private:
    zmq::context_t context_;
    zmq::socket_t inbox_, outbox_, control_;
    std::size_t workgroupSize_, maxBunchSize_;

    SET_LOGGER("I3CLSimClient");

    std::future<void> clientTask_;
};

I3CLSimClient::impl::impl(const std::string &server_address) : context_(1),
    inbox_(context_, ZMQ_PUSH), outbox_(context_, ZMQ_PULL), control_(context_, ZMQ_PUB),
    workgroupSize_(0), maxBunchSize_(0)
{
    inbox_.setsockopt(ZMQ_SNDHWM, 1);
    outbox_.setsockopt(ZMQ_RCVHWM, 1);
    inbox_.bind("inproc://inbox");
    outbox_.bind("inproc://outbox");
    control_.bind("inproc://control");
    
    zmq::socket_t heartbeat(context_, ZMQ_REP);
    heartbeat.bind("inproc://heartbeat");
    clientTask_ = std::async(std::launch::async, ClientWorker, std::ref(context_), server_address, 1);
    
    {
        zmq::message_t address;
        std::list<zmq::message_t> body;
        std::tie(address, body) = read_message(heartbeat);
        {
            zmq::message_t greeting;
            greeting.copy(&I3CLSimServerUtils::greeting);
            heartbeat.send(std::move(greeting));
        }
        i3_assert( body.size() == 1 );
        std::pair<std::size_t,std::size_t> bunchSize = deserialize<std::pair<std::size_t,std::size_t> >(body.front());
        workgroupSize_ = bunchSize.first;
        maxBunchSize_ = bunchSize.second;
    }
    log_trace_stream("connected! granularity "<<workgroupSize_<<" maxSize "<<maxBunchSize_);
}

I3CLSimClient::impl::~impl()
{
    control_.send(zmq::message_t("kbye", 4));
    try {
        if (clientTask_.wait_for(std::chrono::seconds(1)) != std::future_status::ready) {
            log_error_stream("Client task failed to exit");
        }
        // TODO: reap here
    } catch (...) {
        log_error_stream("something went horribly wrong");
    }
}

void I3CLSimClient::impl::EnqueueSteps(I3CLSimStepSeriesConstPtr steps, uint32_t identifier)
{
    log_trace_stream("enqueue steps");
    inbox_.send(serialize(steps), ZMQ_SNDMORE);
    inbox_.send(serialize(identifier), 0);
}

void I3CLSimClient::impl::EnqueueBarrier()
{
    log_trace_stream("enqueue barrier");
    zmq::message_t barrier;
    barrier.copy(&I3CLSimServerUtils::barrier);
    inbox_.send(std::move(barrier));
}

I3CLSimStepToPhotonConverter::ConversionResult_t
I3CLSimClient::impl::GetConversionResultWithBarrierInfo(bool &barrierWasJustReset)
{
    I3CLSimStepToPhotonConverter::ConversionResult_t result;
    std::deque<zmq::message_t> body;
    do {
        zmq::message_t msg;
        outbox_.recv(&msg);
        body.emplace_back(std::move(msg));
    } while (body.back().more());

    // End of message packet is always the identifier
    i3_assert( body.size() >= 2 );
    result.identifier = deserialize<decltype(result.identifier)>(body.back());
    body.pop_back();
    // Start of message packet is always the barrier info
    barrierWasJustReset = deserialize<bool>(body.front());
    body.pop_front();
    
    // Optional: photons
    if (body.size() > 0) {
        result.photons = deserialize<decltype(result.photons)>(body.front());
        body.pop_front();
    }
    // Optional: photon histories
    if (body.size() > 0) {
        result.photonHistories = deserialize<decltype(result.photonHistories)>(body.front());
        body.pop_front();
    }
    i3_assert( body.size() == 0 );
    
    return result;
}

#else // #defined HAS_ZMQ

namespace {

struct no_zmq {
    no_zmq() { throw std::runtime_error("CLSim was built without ZMQ support. I3CLSimClient cannot function."); }
};

}

class I3CLSimServer::impl : no_zmq {
public:
    impl(const std::string &address, const std::vector<I3CLSimStepToPhotonConverterPtr> &converters) {}

    std::map<std::string, double> GetStatistics() const { return std::map<std::string, double>(); }
    std::string GetAddress() const { return ""; }
};

class I3CLSimClient::impl : no_zmq {
public:
    impl(const std::string &server_address) {};

    void EnqueueSteps(I3CLSimStepSeriesConstPtr steps, uint32_t identifier) {}
    void EnqueueBarrier() {}
    I3CLSimStepToPhotonConverter::ConversionResult_t GetConversionResultWithBarrierInfo(bool &barrierWasJustReset) { return I3CLSimStepToPhotonConverter::ConversionResult_t(); }
    std::size_t GetWorkgroupSize() const { return 0; }
    std::size_t GetMaxNumWorkitems() const { return 0; }
};

#endif // #defined HAS_ZMQ

///////////////////////// I3CLSimServer interface /////////////////////////////

I3CLSimServer::I3CLSimServer(
    const std::string &address,
    const std::vector<I3CLSimStepToPhotonConverterPtr> &converters) : 
    impl_(new impl(address, converters))
{}

// Define special functions here, where the implementation type is complete
I3CLSimServer::~I3CLSimServer() = default;
I3CLSimServer::I3CLSimServer(I3CLSimServer &&) = default;
I3CLSimServer& I3CLSimServer::operator=(I3CLSimServer &&) = default;

std::map<std::string, double> I3CLSimServer::GetStatistics() const
{
    return impl_->GetStatistics();
}

std::string I3CLSimServer::GetAddress() const
{
    return impl_->GetAddress();
}

///////////////////////// I3CLSimCleint interface /////////////////////////////

I3CLSimClient::I3CLSimClient(const std::string &server_address)
    : impl_(new impl(server_address))
{}

// Define special functions here, where the implementation type is complete
I3CLSimClient::~I3CLSimClient() = default;
I3CLSimClient::I3CLSimClient(I3CLSimClient &&) = default;
I3CLSimClient& I3CLSimClient::operator=(I3CLSimClient &&) = default;

std::size_t I3CLSimClient::GetWorkgroupSize() const { return impl_->GetWorkgroupSize(); }
std::size_t I3CLSimClient::GetMaxNumWorkitems() const { return impl_->GetMaxNumWorkitems(); }

void I3CLSimClient::EnqueueSteps(I3CLSimStepSeriesConstPtr steps, uint32_t identifier)
{
    impl_->EnqueueSteps(steps, identifier);
}

void I3CLSimClient::EnqueueBarrier()
{
    impl_->EnqueueBarrier();
}

I3CLSimStepToPhotonConverter::ConversionResult_t
I3CLSimClient::GetConversionResultWithBarrierInfo(bool &barrierWasJustReset)
{
    return impl_->GetConversionResultWithBarrierInfo(barrierWasJustReset);
}