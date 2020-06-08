
i3_project(cmc
  PYTHON_DIR python
  DOCS_DIR resources/docs
  )

i3_add_library(cmc 
  private/cmc/I3CascadeMCCommon.cxx
  private/cmc/I3CascadeMCService.cxx
  private/cmc/I3CascadeSplit.cxx
  private/cmc/I3CascadeMuonSplit.cxx
  private/cmc/I3CascadeDevelopment.cxx
  private/cmc/I3CascadeParametrization.cxx
  private/cmc/I3CascadeSimulation.cxx
  private/cmc/I3CascadeSimulationCrossSections.cxx
  private/cmc/I3MetropolisHastings.cxx

  USE_TOOLS python gsl boost

  USE_PROJECTS serialization icetray dataclasses phys-services sim-services
)

i3_add_pybindings(cmc
  private/pybindings/I3CascadeMCService.cxx
  private/pybindings/module.cxx
  USE_TOOLS boost python
  USE_PROJECTS icetray dataclasses cmc
)

set(LIB_${PROJECT_NAME}_TESTS
  private/test/main.cxx
)

if (SPRNG_FOUND)
  # this test requires SPRNG
  LIST(APPEND LIB_${PROJECT_NAME}_TESTS
    private/test/I3CascadeSplitTests.cxx
    private/test/RepeatabilityTests.cxx
  )
endif (SPRNG_FOUND)

i3_test_executable(test
  ${LIB_${PROJECT_NAME}_TESTS}
  USE_TOOLS boost
  USE_PROJECTS serialization dataclasses cmc sim-services dataio icetray phys-services
)

i3_test_scripts(resources/test/test_shape.py)

