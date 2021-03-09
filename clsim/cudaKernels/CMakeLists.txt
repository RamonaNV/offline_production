# /*The MIT License (MIT)

# Copyright (c) 2020, Hendrik Schwanekamp hschwanekamp@nvidia.com, Ramona Hohl rhohl@nvidia.com

# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGSEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
# */
add_library(thiscudalib SHARED
  propagationKernelSource.cu
  I3CLSimCUDADevice.cu
  )

set(_architectures
  35 # K20
  37 # K80
  52 # GTX 980
  61 # P4, GTX 1080
  75 # RTX 2080
)

if (CUDAToolkit_VERSION_MAJOR VERSION_GREATER_EQUAL 11)
  list(APPEND _architectures
    80 # A100
    86 # RTX 3080
  )
elseif ()
  list(APPEND _architectures
    30 # GTX 680
  )
endif()

set(CLSIM_CUDA_ARCHITECTURES ${_architectures} CACHE STRING "CUDA architectures to build for (see https://en.wikipedia.org/wiki/CUDA#GPUs_supported)")

colormsg(GREEN "+-- CUDA architectures: ${CLSIM_CUDA_ARCHITECTURES}")

set_property(TARGET thiscudalib PROPERTY CUDA_ARCHITECTURES ${CLSIM_CUDA_ARCHITECTURES})

# call constexpr __host__ functions from __device__
target_compile_options(thiscudalib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)
target_compile_options(thiscudalib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-use_fast_math>)
# source code information for profiling
target_compile_options(thiscudalib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-lineinfo>)
# forward unkown linker options to host linker to fix problem with older CMake versions (GCC linker flags being passed to nvcc)
target_compile_options(thiscudalib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-forward-unknown-to-host-linker>)


target_include_directories(thiscudalib
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../private
  ${CMAKE_CURRENT_SOURCE_DIR}/../public
  ${CMAKE_CURRENT_SOURCE_DIR}/../private/clsim
  ${CMAKE_CURRENT_SOURCE_DIR}/../public/clsim
  ${CMAKE_CURRENT_SOURCE_DIR}/../../dataclasses/public
  ${CMAKE_CURRENT_SOURCE_DIR}/../../icetray/public
  ${CMAKE_CURRENT_SOURCE_DIR}/../../serialization/public
  ${CMAKE_CURRENT_SOURCE_DIR}/../../phys-services/public
  ${BOOST_INCLUDE_DIR}  # required when building against CVMFS boost using the host CUDA/nvcc
  )
 
# set some options for benchmarking
OPTION(USE_JOBQUEUE "use the load balancing jobqueue" OFF)
IF(USE_JOBQUEUE)
  target_compile_definitions(thiscudalib PUBLIC -DUSE_JOBQUEUE)
ENDIF(USE_JOBQUEUE)

OPTION(BLOCK_RANDOM_NUMBERS "use blocked random numbers" OFF)
IF(BLOCK_RANDOM_NUMBERS)
  target_compile_definitions(thiscudalib PUBLIC -DBLOCK_RANDOM_NUMBERS_SCATTERING -DBLOCK_RANDOM_NUMBERS_PROPAGATION)
ENDIF(BLOCK_RANDOM_NUMBERS)

OPTION(BENCHMARK_OPENCL "run opencl benchmark" OFF)
IF(BENCHMARK_OPENCL)
  target_compile_definitions(thiscudalib PUBLIC -DBENCHMARK_OPENCL)
ENDIF(BENCHMARK_OPENCL)

OPTION(BENCHMARK_OPENCL_SHUFFLED "run opencl shuffled benchmark" OFF)
IF(BENCHMARK_OPENCL_SHUFFLED)
  target_compile_definitions(thiscudalib PUBLIC -DBENCHMARK_OPENCL_SHUFFLED)
ENDIF(BENCHMARK_OPENCL_SHUFFLED)

OPTION(BENCHMARK_OPENCL_SHUFFLED_32 "run opencl shuffled 32 benchmark" OFF)
IF(BENCHMARK_OPENCL_SHUFFLED_32)
  target_compile_definitions(thiscudalib PUBLIC -DBENCHMARK_OPENCL_SHUFFLED_32)
ENDIF(BENCHMARK_OPENCL_SHUFFLED_32)

OPTION(BENCHMARK_SHUFFLED "run shuffled benchmark" OFF)
IF(BENCHMARK_SHUFFLED)
  target_compile_definitions(thiscudalib PUBLIC -DBENCHMARK_SHUFFLED)
ENDIF(BENCHMARK_SHUFFLED)

OPTION(BENCHMARK_SHUFFLED_32 "run shuffled 32 benchmark" OFF)
IF(BENCHMARK_SHUFFLED_32)
  target_compile_definitions(thiscudalib PUBLIC -DBENCHMARK_SHUFFLED_32)
ENDIF(BENCHMARK_SHUFFLED_32)
