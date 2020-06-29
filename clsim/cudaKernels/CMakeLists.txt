
find_package(CUDAToolkit REQUIRED)

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
