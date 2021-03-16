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

find_package(CUDAToolkit REQUIRED)
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS}  -arch=sm_75 -use_fast_math --expt-relaxed-constexpr --keep  --ptxas-options -v  --resource-usage")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS}  -lineinfo") # source code information for profiling


#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake")


find_package(OptiX7 REQUIRED)
include_directories(${OptiX7_INCLUDE_DIRS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/optixSrc)



include(PTXUtilities)
add_library(thiscudalib SHARED
  propagationKernelSource.cu
  optixSrc/rtxFunctions.cpp
  optixSrc/common.cpp
  )

target_link_libraries(thiscudalib
  PUBLIC
  CUDA::cudart
  OptiX7
  )
 

 # add_ptx_targets(thiscudalib optixKernels)
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
  ${CMAKE_CURRENT_SOURCE_DIR}/optixSrc
  )
   
  add_subdirectory(./optixSrc)
