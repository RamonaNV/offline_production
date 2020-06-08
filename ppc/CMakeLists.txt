if(OPENCL_FOUND)
  i3_project(ppc
    PYTHON_DIR python
    DOCS_DIR resources/docs
    )

  colormsg(GREEN "+-- OpenCL SDK found.  Building ppc.")
  i3_add_library(ppc
    private/ppc/i3ppc.cxx
    private/ppc/I3CLSimStepToPhotonConverterPPC.cxx
    USE_TOOLS boost python
    USE_PROJECTS serialization icetray dataclasses simclasses phys-services interfaces clsim
    )
  i3_add_library(xppc
    private/ppc/ocl/ppc.cxx
    USE_TOOLS opencl
    )
  target_link_libraries(ppc xppc)
  target_link_libraries(xppc icetray)
  target_compile_definitions(xppc PRIVATE -DXLIB -DUSE_I3_LOGGING)
  set_target_properties(xppc PROPERTIES
    CXX_EXTENSIONS YES
  )
  i3_add_pybindings(ppc
    private/ppc/pybindings.cxx
    USE_PROJECTS clsim ppc
  )
else()
  message(WARNING "OpenCL was not found.  Not building ppc.")
  find_package(OpenCL) # find non-quiet to let the user know what went wrong.
endif()
