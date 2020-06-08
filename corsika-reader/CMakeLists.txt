i3_project(corsika-reader
	DOCS_DIR resources/docs
	PYTHON_DIR python
)

if (Boost_VERSION VERSION_GREATER 106399)
  MESSAGE(STATUS "    + Boost >= 6.4 : Adding CORSIKA Client/Server")
  ADD_DEFINITIONS(-DUSE_CORSIKA_CLIENT)
  set(PROCESS_FILES 
      private/corsika-reader/CorsikaClient.cxx
      private/corsika-reader/I3CORSIKAService.cxx	
      )
else (Boost_VERSION VERSION_GREATER 106399)
  MESSAGE(STATUS "    - Boost < 6.4 : Skipping CORSIKA Client/Server")
endif (Boost_VERSION VERSION_GREATER 106399)

i3_add_library(corsika-reader
	private/corsika-reader/I3CORSIKAReader.cxx
	private/corsika-reader/CORSIKAResampler.cxx
	private/corsika-reader/I3CORSIKAReaderUtils.cxx
	${PROCESS_FILES}
	USE_TOOLS boost python
	USE_PROJECTS serialization icetray dataio dataclasses simclasses sim-services phys-services 
)

i3_add_pybindings(corsika_reader
      private/pybindings/corsika-reader.cxx
      USE_TOOLS boost python
      USE_PROJECTS corsika-reader
      ) 

i3_executable_script(check resources/scripts/check.py)

#i3_executable(cli private/corsika-reader/cli.cxx USE_TOOLS boost USE_PROJECTS corsika-reader neutrino-generator earthmodel-service MuonGun)

# boost::process wasn't added until 1.64. to build against an earlier version,
# copy or symlink:
# boost/
#   asio
#   process
#   winapi
#   asio.hpp
#   process.hpp
# from the include directory of a newer boost release to private/vendor
#if (Boost_VERSION VERSION_LESS 106400)
#  target_include_directories(corsika-reader PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/private/vendor)
#endif (Boost_VERSION VERSION_LESS 106400)

#target_include_directories(corsika-reader-cli PRIVATE ${CMAKE_SOURCE_DIR}/weighting/public)

i3_test_scripts(
	resources/test/testcorsika.py
)

