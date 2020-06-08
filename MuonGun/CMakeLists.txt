SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wno-inconsistent-missing-override")

i3_project(MuonGun
    PYTHON_DIR python
    DOCS_DIR resources/docs
)

SET (LIB_${PROJECT_NAME}_SOURCEFILES
    # Numerical integration
    private/cubature/cubature.c
	
    # Base classes and function implementations
    private/MuonGun/I3MuonGun.cxx
    private/MuonGun/SplineTable.cxx
    private/MuonGun/Track.cxx
    private/MuonGun/Generator.cxx
    private/MuonGun/WeightCalculator.cxx
    private/MuonGun/SamplingSurface.cxx
    private/MuonGun/Cylinder.cxx
    private/MuonGun/ExtrudedPolygon.cxx
    private/MuonGun/Flux.cxx
    private/MuonGun/RadialDistribution.cxx
    private/MuonGun/EnergyDistribution.cxx
    
    # Event generator implementations
    private/MuonGun/Floodlight.cxx
    private/MuonGun/StaticSurfaceInjector.cxx
    private/MuonGun/EnergyDependentSurfaceInjector.cxx
    private/MuonGun/CORSIKAGenerationProbability.cxx
    private/MuonGun/NaturalRateInjector.cxx
)
SET (LIB_${PROJECT_NAME}_PROJECTS
    serialization icetray dataclasses phys-services simclasses photospline tableio
)
SET (LIB_${PROJECT_NAME}_TOOLS
    boost python gsl
)

SET (LIB_${PROJECT_NAME}_PYSOURCES 
    private/pybindings/I3MuonGun.cxx
    private/pybindings/Track.cxx
    private/pybindings/Injectors.cxx
    private/pybindings/Generator.cxx
    private/pybindings/Surface.cxx
    private/pybindings/Flux.cxx
    private/pybindings/RadialDistribution.cxx
    private/pybindings/EnergyDistribution.cxx
    private/pybindings/WeightCalculator.cxx
    private/pybindings/CORSIKAGenerationProbability.cxx
    private/pybindings/module.cxx
)

SET(BUILD_MUONITRON YES)
IF (NOT IS_DIRECTORY ${CMAKE_SOURCE_DIR}/PROPOSAL)
    COLORMSG (CYAN "+-- PROPOSAL not found, disabling Muonitron table-maker")
    SET(BUILD_MUONITRON NO)
ENDIF()
IF (${NUMPY_FOUND})
    LIST(APPEND LIB_${PROJECT_NAME}_TOOLS numpy)
ELSE (${NUMPY_FOUND})
    COLORMSG (CYAN "+-- Numpy not found, disabling Muonitron table-maker")
    SET(BUILD_MUONITRON NO)
ENDIF (${NUMPY_FOUND})

IF (${BUILD_MUONITRON})
    COLORMSG (GREEN "+-- Muonitron table-maker")
    LIST (APPEND LIB_${PROJECT_NAME}_SOURCEFILES
        # Support bits for making muon-flux parametrizations
        private/MuonGun/Muonitron.cxx
        private/MuonGun/CompactTrack.cxx
        private/MuonGun/TrackBinner.cxx
        private/MuonGun/MuonPropagator.cxx
        private/MuonGun/histogram.cxx
    )
    LIST (APPEND LIB_${PROJECT_NAME}_PYSOURCES
        private/pybindings/CompactTrack.cxx
        private/pybindings/TrackBinner.cxx
        private/pybindings/MuonPropagator.cxx
        private/pybindings/histogram.cxx
    )
    LIST (APPEND LIB_${PROJECT_NAME}_PROJECTS
        PROPOSAL
    )
    ADD_DEFINITIONS(-DUSE_PROPOSAL)
ENDIF (${BUILD_MUONITRON})

i3_add_library(MuonGun
    ${LIB_${PROJECT_NAME}_SOURCEFILES}
    USE_TOOLS ${LIB_${PROJECT_NAME}_TOOLS}
    USE_PROJECTS ${LIB_${PROJECT_NAME}_PROJECTS}
)

IF (BOOSTNUMPY_FOUND)
    LIST (APPEND LIB_${PROJECT_NAME}_TOOLS boostnumpy)
    ADD_DEFINITIONS(-DUSE_NUMPY)
    COLORMSG (CYAN "+-- boost::numpy")
ENDIF (BOOSTNUMPY_FOUND)

i3_add_pybindings(MuonGun
  ${LIB_${PROJECT_NAME}_PYSOURCES}
  USE_TOOLS ${LIB_${PROJECT_NAME}_TOOLS}
  USE_PROJECTS MuonGun
  )

i3_test_scripts(
  resources/test/test_StaticSurfaceInjector.py
  resources/test/test_NaturalRateInjector.py
  resources/test/extruded_polygon.py
)

i3_test_executable(test 
  private/test/main.cxx
  private/test/SplineTable.cxx
  private/test/Surface.cxx
  private/test/EnergyDistribution.cxx
  private/test/Generator.cxx
  private/test/Integration.cxx
  private/test/EnsembleSampler.cxx
  USE_PROJECTS MuonGun icetray dataclasses phys-services
)

