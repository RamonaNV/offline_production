#
#  $Id: CMakeLists.txt 169307 2016-06-01 18:27:08Z cweaver $
#
i3_project(interfaces
  DOCS_DIR resources/docs
)

i3_add_library(interfaces
  private/interfaces/*.cxx
  USE_PROJECTS serialization icetray
  USE_TOOLS boost
  )

i3_add_pybindings(interfaces
  private/pybindings/module.cxx
  private/pybindings/I3GeometryService.cxx
  private/pybindings/I3CalibrationService.cxx
  private/pybindings/I3DetectorStatusService.cxx
  private/pybindings/I3EventService.cxx
  private/pybindings/I3MetaService.cxx
  private/pybindings/I3OMKey2MBID.cxx
  USE_TOOLS boost python
  USE_PROJECTS interfaces dataclasses
  )


