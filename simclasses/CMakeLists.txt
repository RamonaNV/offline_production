
i3_project(simclasses
  PYTHON_DIR python
  DOCS_DIR resources/docs)

SET(LIB_${PROJECT_NAME}_SOURCEFILES 
  private/simclasses/CorsikaLongStep.cxx
  private/simclasses/I3MMCTrack.cxx
  private/simclasses/I3MCNKGInfo.cxx
  private/simclasses/I3CorsikaShowerInfo.cxx
  private/simclasses/I3MCPulse.cxx
  private/simclasses/I3MCPE.cxx
  private/simclasses/I3WimpParams.cxx
  private/simclasses/I3ParticleIDMap.cxx
  private/simclasses/I3Photon.cxx
  private/simclasses/I3CompressedPhoton.cxx
  private/simclasses/I3CylinderMap.cxx
  private/simclasses/I3ExtraGeometryItem.cxx   
  private/simclasses/I3ExtraGeometryItemUnion.cxx   
  private/simclasses/I3ExtraGeometryItemMove.cxx   
  private/simclasses/I3ExtraGeometryItemCylinder.cxx
  private/simclasses/I3NuGenInfo.cxx  
  private/simclasses/I3CorsikaInfo.cxx  
  private/simclasses/I3ShowerBias.cxx
  private/simclasses/I3CorsikaWeight.cxx
  
  # tableio converters
  private/simclasses/converter/I3PhotonConverter.cxx
  private/simclasses/converter/I3MCPEConverterWithIDs.cxx
  
  # the following are deprecated, but need to be 
  # kept indefinitely so we can still read old data
  private/deprecated/I3MCPMTResponse.cxx
  private/deprecated/I3GaussianPMTPulse.cxx
  private/deprecated/I3MCTWRParams.cxx
  private/deprecated/MMCWeight.cxx
  )


SET(LIB_${PROJECT_NAME}_PROJECTS serialization icetray dataclasses phys-services)
SET(LIB_${PROJECT_NAME}_TOOLS boost python)

if(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/tableio)
	LIST(APPEND LIB_${PROJECT_NAME}_PROJECTS tableio)
	LIST(APPEND LIB_${PROJECT_NAME}_SOURCEFILES private/simclasses/converter/*.cxx)
	ADD_DEFINITIONS(-DUSE_TABLEIO)
	MESSAGE(STATUS "     + tableio support")
endif(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/tableio)

i3_add_library(${PROJECT_NAME} ${LIB_${PROJECT_NAME}_SOURCEFILES}
	USE_TOOLS ${LIB_${PROJECT_NAME}_TOOLS}
	USE_PROJECTS ${LIB_${PROJECT_NAME}_PROJECTS}
)

i3_add_pybindings(simclasses
  private/pybindings/CorsikaLongStep.cxx
  private/pybindings/I3Converters.cxx
  private/pybindings/I3CorsikaShowerInfo.cxx
  private/pybindings/I3MMCTrack.cxx
  private/pybindings/I3MCPulse.cxx
  private/pybindings/I3ParticleIDMap.cxx
  private/pybindings/I3MCPE.cxx
  private/pybindings/I3Photon.cxx
  private/pybindings/I3CompressedPhoton.cxx
  private/pybindings/I3WimpParams.cxx  
  private/pybindings/module.cxx
  private/pybindings/I3CylinderMap.cxx
  private/pybindings/I3ExtraGeometryItem.cxx
  private/pybindings/const_ptr_helpers.h  
  private/pybindings/python_gil_holder.h 
  private/pybindings/I3NuGenInfo.cxx 
  private/pybindings/I3CorsikaInfo.cxx 
  private/pybindings/I3CorsikaWeight.cxx
  private/pybindings/I3ShowerBias.cxx  
  USE_TOOLS boost python
  USE_PROJECTS simclasses icetray tableio
  )

i3_test_scripts(resources/test/*.py)
