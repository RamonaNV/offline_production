
#include "ppc/I3CLSimStepToPhotonConverterPPC.h"

#include <icetray/load_project.h>
namespace bp = boost::python;

static void register_I3CLSimStepToPhotonConverterPPC()
{
    bp::class_<I3CLSimStepToPhotonConverterPPC,
               boost::shared_ptr<I3CLSimStepToPhotonConverterPPC>,
               bp::bases<I3CLSimStepToPhotonConverter>,
               boost::noncopyable >("I3CLSimStepToPhotonConverterPPC", bp::no_init)
        .add_static_property("instance", &I3CLSimStepToPhotonConverterPPC::GetInstance)
    ;
}

BOOST_PYTHON_MODULE(ppc)
{
    load_project("ppc", false);
    
    bp::import("icecube.clsim");
    
    register_I3CLSimStepToPhotonConverterPPC();
}
