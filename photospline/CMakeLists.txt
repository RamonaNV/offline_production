#
# $Id: CMakeLists.txt 64301 2010-05-25 16:49:40Z jvansanten $
#
set(OLD_DOXYGEN_FOUND ${DOXYGEN_FOUND})
unset(DOXYGEN_FOUND)
i3_project(photospline
	PYTHON_DIR python
        DOCS_DIR resources/docs
)

if(BLAS_FOUND)
	set(PHOTOSPLINE_EXTRA_TOOLS blas)
endif(BLAS_FOUND)

i3_add_library(photospline
	private/lib/bspline.c
	private/lib/bspline_multi.c
	private/lib/convolve.c
	private/lib/fitstable.c
	private/lib/splinepdf.c
	private/lib/I3SplineTable.cxx
	USE_TOOLS cfitsio gsl python ${PHOTOSPLINE_EXTRA_TOOLS}
)

if(NUMPY_FOUND)
i3_add_pybindings(photospline
	private/pybindings/I3SplineTable.cxx
	private/pybindings/module.cxx
	USE_TOOLS boost python numpy
	USE_PROJECTS icetray photospline
)
endif(NUMPY_FOUND)

i3_executable(evalsplinefits
	private/util/evalsplinefits.c
	USE_TOOLS cfitsio
	USE_PROJECTS photospline)

i3_executable(dumpsplinefits
	private/util/dumpsplinefits.c
	USE_TOOLS cfitsio
	USE_PROJECTS photospline)

SET_TARGET_PROPERTIES(photospline-evalsplinefits
        PROPERTIES
        COMPILE_FLAGS "-std=c99"
)

if(NOT DEFINED PHOTOSPLINE_CFLAGS)
	# all the goodies of -ffast-math, minus -funsafe-math-optimizations
	SET(_ffast_math "-fno-math-errno -fno-trapping-math -ffinite-math-only")
	# only use fno-signaling-nans with gcc
	if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
		SET(_ffast_math "${_ffast_math} -fno-signaling-nans")
	endif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
	SET(PHOTOSPLINE_CFLAGS "-O3 ${_ffast_math} -std=c99")
	# fpeel-loops is only valid for gcc, not clang or pgcc
	if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
		SET(PHOTOSPLINE_CFLAGS "${PHOTOSPLINE_CFLAGS} -fpeel-loops")
	endif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
	# add flags to enable the vector units present on the host processor
	# we'd love to use -march=native, but that's only available on
	# GCC >= 4.2 and only on some architectures.
	if(CMAKE_SYSTEM_PROCESSOR MATCHES "^i.86$")
		SET(PHOTOSPLINE_CFLAGS "${PHOTOSPLINE_CFLAGS} -msse2 -mfpmath=sse")
	elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(powerpc|ppc)")
		SET(PHOTOSPLINE_CFLAGS "${PHOTOSPLINE_CFLAGS} -maltivec")
	elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^sparc")
		SET(PHOTOSPLINE_CFLAGS "${PHOTOSPLINE_CFLAGS} -mvis")
	endif(CMAKE_SYSTEM_PROCESSOR MATCHES "^i.86$")
	# Note: SSE2 is a requirement of amd64, so we don't worry about that
	SET(PHOTOSPLINE_CFLAGS ${PHOTOSPLINE_CFLAGS} CACHE STRING "Compiler flags specific to photospline")
endif(NOT DEFINED PHOTOSPLINE_CFLAGS)

# We wants it fast, precious.
SET_SOURCE_FILES_PROPERTIES(
	private/lib/bspline.c
	private/lib/bspline_multi.c
	private/lib/convolve.c
	private/lib/fitstable.c
	private/lib/splinepdf.c
	PROPERTIES
	COMPILE_FLAGS ${PHOTOSPLINE_CFLAGS}
)

i3_test_executable(test 
	private/test/*.cxx
	USE_TOOLS boost
	USE_PROJECTS photospline icetray
)

execute_process(COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/resources/test/find_pyfits.py" RESULT_VARIABLE PYFITS_FOUND ERROR_QUIET)
# we need a NOT here because PYFITS_FOUND is 0 if the process succeeded
if (NOT PYFITS_FOUND)
	list(APPEND _test_scripts resources/test/deriv2.py)
endif (NOT PYFITS_FOUND)

include(CheckIncludeFile)

if(NUMPY_FOUND AND BLAS_FOUND AND LAPACK_FOUND AND SUITESPARSE_FOUND)
	i3_add_library(spglam
		private/cfitter/cholesky_solve.c private/cfitter/glam.c
		private/cfitter/nnls.c private/cfitter/pyglam.c
		private/cfitter/splineutil.c
		USE_TOOLS blas lapack python numpy suitesparse
		USE_PROJECTS photospline
		INSTALL_DESTINATION lib/icecube/photospline
		NOT_INSPECTABLE
		MODULE
	)
	if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
		target_link_libraries(spglam rt)
	endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	# SPQR pulls in libc++ dependencies to an otherwise pure-C target.
	# Explicitly link it on El Capitan and greater
	if(APPLE AND CMAKE_SYSTEM_VERSION VERSION_GREATER 14)
		target_link_libraries(spglam c++)
	endif(APPLE AND CMAKE_SYSTEM_VERSION VERSION_GREATER 14)
	add_custom_command(TARGET spglam
		PRE_LINK
		COMMAND mkdir -p ${CMAKE_BINARY_DIR}/lib/icecube/photospline
	)
	set_target_properties(spglam PROPERTIES
		PREFIX ""
		OUTPUT_NAME spglam 
		LIBRARY_OUTPUT_DIRECTORY
		${LIBRARY_OUTPUT_PATH}/icecube/photospline
	)
	add_custom_command(TARGET spglam
		PRE_LINK
		COMMAND mkdir -p ${CMAKE_BINARY_DIR}/lib/icecube/photospline
	)
	execute_process(COMMAND ${PYTHON_EXECUTABLE} "-c \"import astropy\"" RESULT_VARIABLE ASTROPY_FOUND ERROR_QUIET)
	# we need a NOT here because PYFITS_FOUND is 0 if the process succeeded
	if (NOT ASTROPY_FOUND)
		list(APPEND _test_scripts resources/test/test_cfitter.py)
	endif()
endif(NUMPY_FOUND AND BLAS_FOUND AND LAPACK_FOUND AND SUITESPARSE_FOUND)

i3_test_scripts(${_test_scripts})

set(DOXYGEN_FOUND ${OLD_DOXYGEN_FOUND})
