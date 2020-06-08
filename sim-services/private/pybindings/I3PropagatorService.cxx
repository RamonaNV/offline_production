#include <sim-services/I3PropagatorService.h>
#include <icetray/python/gil_holder.hpp>
#include <icetray/python/std_map_indexing_suite.hpp>

namespace bp = boost::python;

struct I3PropagatorServiceWrapper : I3PropagatorService, bp::wrapper<I3PropagatorService>
{
    // pure virtual
    virtual std::vector<I3Particle> Propagate(I3Particle& p, DiagnosticMapPtr frame, I3FramePtr) {
        // python needs to be able to change this, so provide it with a pointer
        I3ParticlePtr p_(new I3Particle(p));

        boost::python::detail::gil_holder gil;
        std::vector<I3Particle> children = 
            this->get_override("Propagate")(p_, frame);

        p = *p_;
        return children;
    }

    // pure virtual
    virtual void SetRandomNumberGenerator(I3RandomServicePtr random)
    {
        boost::python::detail::gil_holder gil;
        this->get_override("SetRandomNumberGenerator")(random);
    }

};

void register_I3PropagatorService()
{
    {
        bp::scope I3PropagatorService_scope = 
        bp::class_<I3PropagatorServiceWrapper, boost::shared_ptr<I3PropagatorServiceWrapper>, boost::noncopyable>("I3PropagatorService")
        .def("Propagate", bp::pure_virtual(&I3PropagatorService::Propagate), (bp::args("particle"), bp::args("diagnostics")=boost::shared_ptr<I3PropagatorService::DiagnosticMap>(),bp::args("frame")=boost::shared_ptr<I3Frame>()))
        .def("SetRandomNumberGenerator", bp::pure_virtual(&I3PropagatorService::SetRandomNumberGenerator))
        ;
        
        bp::class_<I3PropagatorService::DiagnosticMap, boost::shared_ptr<I3PropagatorService::DiagnosticMap> >("DiagnosticMap", bp::no_init)
        ;
    }
    
    bp::implicitly_convertible<boost::shared_ptr<I3PropagatorServiceWrapper>, boost::shared_ptr<const I3PropagatorService> >();
    bp::implicitly_convertible<boost::shared_ptr<I3PropagatorServiceWrapper>, boost::shared_ptr<I3PropagatorService> >();
    bp::implicitly_convertible<boost::shared_ptr<I3PropagatorServiceWrapper>, boost::shared_ptr<const I3PropagatorServiceWrapper> >();

    bp::class_<I3ParticleTypePropagatorServiceMap, I3ParticleTypePropagatorServiceMapPtr>("I3ParticleTypePropagatorServiceMap")
    .def(bp::std_map_indexing_suite<I3ParticleTypePropagatorServiceMap>())
    ;
}
