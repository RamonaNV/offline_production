
#include "clsim/dom/I3CLSimPhotonToMCPEConverter.h"
#include "clsim/dom/I3PhotonToMCPEConverter.h"
#include "simclasses/I3CompressedPhoton.h"

namespace bp = boost::python;

bp::object Convert_Photon(const I3CLSimPhotonToMCPEConverter &self, const ModuleKey &key, const I3CompressedPhoton &photon)
{
    auto hit = self.Convert(key, photon);
    if (hit) {
        return bp::make_tuple(std::get<0>(*hit), std::get<1>(*hit));
    } else {
        return bp::object();
    }
}

void register_I3CLSimPhotonToMCPEConverter()
{
    bp::class_<I3CLSimPhotonToMCPEConverter, boost::shared_ptr<I3CLSimPhotonToMCPEConverter>, boost::noncopyable>("I3CLSimPhotonToMCPEConverter", bp::no_init)
        .def("Convert", Convert_Photon)
    ;
}

void register_I3CLSimPhotonToMCPEConverterForDOMs()
{
    bp::class_<I3CLSimPhotonToMCPEConverterForDOMs, boost::shared_ptr<I3CLSimPhotonToMCPEConverterForDOMs>,
        bp::bases<I3CLSimPhotonToMCPEConverter>, boost::noncopyable>
        ("I3CLSimPhotonToMCPEConverterForDOMs",
          bp::init<I3RandomServicePtr,boost::shared_ptr<const std::map<OMKey, I3CLSimFunctionConstPtr>>,I3CLSimFunctionConstPtr>(
                   bp::args("randomService", "wavelengthAcceptance", "angularAcceptance")))
    ;
}