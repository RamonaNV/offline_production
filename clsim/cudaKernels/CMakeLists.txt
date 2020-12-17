find_package(CUDAToolkit REQUIRED)
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_75,code=sm_75 -gencode arch=compute_80,code=sm_80 --keep  -use_fast_math --expt-relaxed-constexpr --ptxas-options -v --resource-usage")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS}  -lineinfo") # source code information for profiling

add_library(thiscudalib SHARED
  propagationKernelSource.cu
  )

target_link_libraries(thiscudalib
  PUBLIC
  CUDA::cudart
  )

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