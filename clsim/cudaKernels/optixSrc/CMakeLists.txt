 
add_ptx_targets(thiscudalib optixKernels)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}
${CMAKE_CURRENT_SOURCE_DIR}/..
${CMAKE_CURRENT_SOURCE_DIR}/../../private
${CMAKE_CURRENT_SOURCE_DIR}/../../public
${CMAKE_CURRENT_SOURCE_DIR}/../../private/clsim
${CMAKE_CURRENT_SOURCE_DIR}/../../public/clsim
${CMAKE_CURRENT_SOURCE_DIR}/../../../dataclasses/public
${CMAKE_CURRENT_SOURCE_DIR}/../../../icetray/public
${CMAKE_CURRENT_SOURCE_DIR}/../../../serialization/public
${CMAKE_CURRENT_SOURCE_DIR}/../../../phys-services/public
)


