option(GEANT4_DATA "Include Geant4 data in tarball" OFF)

find_library(LZMA_LIBRARIES lzma
      ENV LD_LIBRARY_PATH
)

i3_project(simprod-scripts
  PYTHON_DEST icecube/simprod
  PYTHON_DIR python 
  DOCS_DIR resources/docs
)

if (SPRNG_FOUND)
  set(TEST_SCRIPT_LIST
    resources/tests/segments/*.py
    resources/tests/generators/*.py
  )
endif(SPRNG_FOUND)    

if(NOT ${COVERAGE})
  i3_test_scripts(${TEST_SCRIPT_LIST})
endif()

if(NOT EXISTS resources/testdata)
    file(MAKE_DIRECTORY resources/testdata)
endif(NOT EXISTS resources/testdata)

add_custom_target(production-tarball
  COMMENT "Building production tarball of ${CMAKE_INSTALL_PREFIX}"
)

foreach(arch_lib 
	${LIBARCHIVE_LIBRARIES} 
	${LAPACK_LIBRARIES} 
	${LZMA_LIBRARIES} 
	${HDF5_LIBRARIES} 
	${LHAPDF_LIB}
	${PYTHIA_SO_LIB}
	${GEANT4_LIBRARIES}
   )
   if( NOT ${arch_lib}  STREQUAL "-lgfortran" )
      add_custom_command(TARGET production-tarball
        PRE_BUILD
        COMMAND mkdir -p  ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}/lib/tools
        COMMAND ${CMAKE_SOURCE_DIR}/cmake/install_shlib.py ${arch_lib} ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}/lib/tools
        COMMENT "Install additional libraries for production: ${arch_lib}"
        )
   endif()
endforeach(arch_lib ${LIBARCHIVE_LIBRARIES})


if (GEANT4_FOUND)
  # determine GEANT4 library directory
  foreach(g4lib ${GEANT4_LIBRARIES})
     get_filename_component(GEANT4_LIBDIR ${g4lib} PATH)
     break()
  endforeach(g4lib ${GEANT4_LIBRARIES})

  # add GEANT4 library directory
  add_custom_command(TARGET production-tarball
     PRE_BUILD
     COMMAND cp -r ${GEANT4_LIBDIR} ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}/lib/
     COMMENT "Installing Geant4 libs"
     )
  add_custom_command(TARGET production-tarball
     PRE_BUILD
     COMMAND mkdir -p ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}/bin
     COMMAND cp ${GEANT4_SH} ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}/bin/
     COMMENT "Installing geant4.sh "
     )
endif (GEANT4_FOUND)

if (GEANT4_DATA)
     add_custom_command(TARGET production-tarball
     PRE_BUILD
     COMMAND mkdir -p ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}/share
     COMMAND cp -r ${I3_PORTS}/share/geant4 ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}/share/
     COMMENT "Installing Geant4 data"
)
endif (GEANT4_DATA)

if (BUILD_GENIE-ICETRAY)
add_custom_command(TARGET production-tarball
     PRE_BUILD
     COMMAND mkdir -p ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}/share/genie
     COMMAND cp -r ${GENIE_PATH}/* ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}/share/genie/
     COMMENT "Installing GENIE data"
)
endif (BUILD_GENIE-ICETRAY)

add_dependencies(production-tarball tarball-install)
add_dependencies(tarball-install-tools production-tarball)
