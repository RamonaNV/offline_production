i3_project(sim-services
  PYTHON_DIR python
  PYTHON_DEST icecube/sim_services
  DOCS_DIR resources/docs
)

# define the lists used to build the basic library
SET(LIB_${PROJECT_NAME}_SOURCEFILES
  private/sim-services/I3SumGenerator.cxx
  private/sim-services/I3GeneratorService.cxx
  private/sim-services/I3CosmicEventGenerator.cxx
  private/sim-services/I3MCTreeHybridSimulationSplitter.cxx
  private/sim-services/I3PropagatorService.cxx
  private/sim-services/I3InIceCORSIKATrimmer.cxx
  private/sim-services/I3MCEventHeaderGenerator.cxx
  private/sim-services/I3PropagatorModule.cxx
  private/sim-services/I3CombineMCPE.cxx
  private/sim-services/I3DownsampleMCPE.cxx
  private/sim-services/I3MCPEtoI3MCHitConverter.cxx
  private/sim-services/I3RemoveLargeDT.cxx
  private/sim-services/I3ModifyEventID.cxx # deprecated
  private/sim-services/I3MCPEConverters.cxx
  private/sim-services/MCPEMCPulseTools.cxx
  private/sim-services/I3SimConstants.cxx
  private/sim-services/I3PrimaryPulseMapper.cxx
)
SET(LIB_${PROJECT_NAME}_TOOLS python boost)

i3_add_library(sim-services
  ${LIB_${PROJECT_NAME}_SOURCEFILES}

  USE_TOOLS ${LIB_${PROJECT_NAME}_TOOLS}

  USE_PROJECTS serialization icetray phys-services dataclasses simclasses
)

i3_add_pybindings(sim_services
  private/pybindings/I3CosmicEventGenerator.cxx
  private/pybindings/I3GeneratorService.cxx
  private/pybindings/I3SumGenerator.cxx
  private/pybindings/I3PropagatorService.cxx
  private/pybindings/I3SimConstants.cxx
  private/pybindings/module.cxx
  USE_TOOLS ${LIB_${PROJECT_NAME}_TOOLS}
  USE_PROJECTS sim-services
)

i3_test_scripts(resources/tests/propagator_state_storage.py
  resources/tests/test_combine_PEs.py
  resources/tests/test_down_sample_mcpe.py
  resources/tests/test_mcpe_to_mchit_converter.py
  resources/tests/test_inice_corsika_trimmer.py
  resources/tests/remove_large_dt_pes.py
)

i3_test_executable(test 
  private/test/main.cxx
  private/test/SortMCHitsTest.cxx
  private/test/I3PrimaryPulseMapper.cxx
  USE_PROJECTS icetray simclasses sim-services
  USE_TOOLS boost python
  )
