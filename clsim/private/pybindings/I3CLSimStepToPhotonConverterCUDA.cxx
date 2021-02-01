#include <cuda/I3CLSimStepToPhotonConverterCUDA.h>

namespace bp = boost::python;

void register_I3CLSimStepToPhotonConverterCUDA()
{
    // I3CLSimStepToPhotonConverterCUDA
    {

        
        bp::class_<
        I3CLSimStepToPhotonConverterCUDA, 
        boost::shared_ptr<I3CLSimStepToPhotonConverterCUDA>, 
        bp::bases<I3CLSimStepToPhotonConverter>,
        boost::noncopyable
        >
        (
         "I3CLSimStepToPhotonConverterCUDA",
         bp::init<
         I3RandomServicePtr,bool
         >(
           (
            bp::arg("RandomService"),
            bp::arg("UseNativeMath")=I3CLSimStepToPhotonConverterCUDA::default_useNativeMath
           )
          )
        )
        .def("Compile", &I3CLSimStepToPhotonConverterCUDA::Compile)
        .def("GetFullSource", &I3CLSimStepToPhotonConverterCUDA::GetFullSource)
        
        .def("GetGeometrySource", &I3CLSimStepToPhotonConverterCUDA::GetGeometrySource)
        .def("GetCollisionDetectionSource", &I3CLSimStepToPhotonConverterCUDA::GetCollisionDetectionSource)

        .def("SetDevice", &I3CLSimStepToPhotonConverterCUDA::SetDevice)
        .def("GetMaxWorkgroupSize", &I3CLSimStepToPhotonConverterCUDA::GetMaxWorkgroupSize)
        .add_property("maxWorkgroupSize", &I3CLSimStepToPhotonConverterCUDA::GetMaxWorkgroupSize)

        .def("GetWorkgroupSize", &I3CLSimStepToPhotonConverterCUDA::GetWorkgroupSize)
        .def("SetWorkgroupSize", &I3CLSimStepToPhotonConverterCUDA::SetWorkgroupSize)
        .def("GetMaxNumWorkitems", &I3CLSimStepToPhotonConverterCUDA::GetMaxNumWorkitems)
        .def("SetMaxNumWorkitems", &I3CLSimStepToPhotonConverterCUDA::SetMaxNumWorkitems)

        .def("SetEnableDoubleBuffering", &I3CLSimStepToPhotonConverterCUDA::SetEnableDoubleBuffering)
        .def("GetEnableDoubleBuffering", &I3CLSimStepToPhotonConverterCUDA::GetEnableDoubleBuffering)

        .def("SetDoublePrecision", &I3CLSimStepToPhotonConverterCUDA::SetDoublePrecision)
        .def("GetDoublePrecision", &I3CLSimStepToPhotonConverterCUDA::GetDoublePrecision)

        .def("SetStopDetectedPhotons", &I3CLSimStepToPhotonConverterCUDA::SetStopDetectedPhotons)
        .def("GetStopDetectedPhotons", &I3CLSimStepToPhotonConverterCUDA::GetStopDetectedPhotons)

        .def("SetSaveAllPhotons", &I3CLSimStepToPhotonConverterCUDA::SetSaveAllPhotons)
        .def("GetSaveAllPhotons", &I3CLSimStepToPhotonConverterCUDA::GetSaveAllPhotons)

        .def("SetSaveAllPhotonsPrescale", &I3CLSimStepToPhotonConverterCUDA::SetSaveAllPhotonsPrescale)
        .def("GetSaveAllPhotonsPrescale", &I3CLSimStepToPhotonConverterCUDA::GetSaveAllPhotonsPrescale)

        .def("SetPhotonHistoryEntries", &I3CLSimStepToPhotonConverterCUDA::SetPhotonHistoryEntries)
        .def("GetPhotonHistoryEntries", &I3CLSimStepToPhotonConverterCUDA::GetPhotonHistoryEntries)

        .def("SetFixedNumberOfAbsorptionLengths", &I3CLSimStepToPhotonConverterCUDA::SetFixedNumberOfAbsorptionLengths)
        .def("GetFixedNumberOfAbsorptionLengths", &I3CLSimStepToPhotonConverterCUDA::GetFixedNumberOfAbsorptionLengths)

        .def("SetDOMPancakeFactor", &I3CLSimStepToPhotonConverterCUDA::SetDOMPancakeFactor)
        .def("GetDOMPancakeFactor", &I3CLSimStepToPhotonConverterCUDA::GetDOMPancakeFactor)


        .add_property("workgroupSize", &I3CLSimStepToPhotonConverterCUDA::GetWorkgroupSize, &I3CLSimStepToPhotonConverterCUDA::SetWorkgroupSize)
        .add_property("maxNumWorkitems", &I3CLSimStepToPhotonConverterCUDA::GetMaxNumWorkitems, &I3CLSimStepToPhotonConverterCUDA::SetMaxNumWorkitems)
        .add_property("enableDoubleBuffering", &I3CLSimStepToPhotonConverterCUDA::GetEnableDoubleBuffering, &I3CLSimStepToPhotonConverterCUDA::SetEnableDoubleBuffering)
        .add_property("doublePrecision", &I3CLSimStepToPhotonConverterCUDA::GetDoublePrecision, &I3CLSimStepToPhotonConverterCUDA::SetDoublePrecision)
        .add_property("stopDetectedPhotons", &I3CLSimStepToPhotonConverterCUDA::GetStopDetectedPhotons, &I3CLSimStepToPhotonConverterCUDA::SetStopDetectedPhotons)
        .add_property("saveAllPhotons", &I3CLSimStepToPhotonConverterCUDA::GetSaveAllPhotons, &I3CLSimStepToPhotonConverterCUDA::SetSaveAllPhotons)
        .add_property("saveAllPhotonsPrescale", &I3CLSimStepToPhotonConverterCUDA::GetSaveAllPhotonsPrescale, &I3CLSimStepToPhotonConverterCUDA::SetSaveAllPhotonsPrescale)
        .add_property("photonHistoryEntries", &I3CLSimStepToPhotonConverterCUDA::GetPhotonHistoryEntries, &I3CLSimStepToPhotonConverterCUDA::SetPhotonHistoryEntries)
        .add_property("fixedNumberOfAbsorptionLengths", &I3CLSimStepToPhotonConverterCUDA::GetFixedNumberOfAbsorptionLengths, &I3CLSimStepToPhotonConverterCUDA::SetFixedNumberOfAbsorptionLengths)
        .add_property("DOMPancakeFactor", &I3CLSimStepToPhotonConverterCUDA::GetDOMPancakeFactor, &I3CLSimStepToPhotonConverterCUDA::SetDOMPancakeFactor)
        ;
    }

    bp::implicitly_convertible<boost::shared_ptr<I3CLSimStepToPhotonConverterCUDA>, boost::shared_ptr<const I3CLSimStepToPhotonConverterCUDA> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimStepToPhotonConverterCUDA>, boost::shared_ptr<I3CLSimStepToPhotonConverter> >();
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimStepToPhotonConverterCUDA>, boost::shared_ptr<const I3CLSimStepToPhotonConverter> >();
}
