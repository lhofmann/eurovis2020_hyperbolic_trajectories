set(module         HyperbolicTrajectories)
set(sources        vtkApproximateDHT.cpp)
set(headers        vtkApproximateDHT.h)
set(module_depends vtkCommonExecutionModel)

if (_legacy_build_system)  
  add_library(${module} ${sources})
  target_link_libraries(${module} PRIVATE ${VTK_LIBRARIES})
  # store header file names in a custom target property (used for generating ParaView plugin)
  set_target_properties(${module} PROPERTIES HEADERS ${headers})
  generate_export_header(${module} EXPORT_FILE_NAME "${module}Module.h")
  target_include_directories(${module} PUBLIC ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
  install(TARGETS ${module} DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

  if (VTK_WRAP_PYTHON)
    if (${VTK_MAJOR_VERSION} GREATER 7)
      # setup VTK module dependencies
      set(MODULE_HIERARCHY_NAME ${module}Hierarchy)
      set(${module}_WRAP_DEPENDS ${module_depends})
      include(vtkWrapHierarchy)
      vtk_wrap_hierarchy(${module} ${CMAKE_CURRENT_BINARY_DIR} ${headers})
      set(KIT_HIERARCHY_FILE ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_HIERARCHY_NAME}.txt)
      set(LIB_HIERARCHY_STAMP ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_HIERARCHY_NAME}.stamp.txt)
    endif ()

    add_library(${module}Python MODULE ${module}PythonInit.cxx)
    set_target_properties(${module}Python PROPERTIES PREFIX "")

    include(vtkWrapPython)
    vtk_wrap_python3(${module}Python ${module}Python_SOURCES "${headers}")

    add_library(${module}PythonD STATIC ${${module}Python_SOURCES})
    set_target_properties(${module}PythonD PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
    target_link_libraries(${module}PythonD ${module})
    target_link_libraries(${module}PythonD
      ${VTK_LIBRARIES}
      vtkWrappingPythonCore
      ${VTK_PYTHON_LIBRARIES})

    # link against PythonD modules that we depend on
    set(_pythond_libs)
    foreach(vtk_module ${module_depends})
      set(_pythond_libs ${_pythond_libs} ${vtk_module}PythonD)
    endforeach()
    target_link_libraries(${module}Python ${module}PythonD ${_pythond_libs})

    set_target_properties(${module}Python PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${_python_module_path}/${_python_module_name})    
    set_target_properties(${module}Python PROPERTIES INSTALL_RPATH "${_python_module_rpath}")
    install(TARGETS ${module}Python DESTINATION ${CMAKE_INSTALL_PREFIX}/${_python_module_path}/${_python_module_name})

    # generate Module.py which imports ModulePython.so
    set(module_py "from __future__ import absolute_import\nfrom .${module}Python import *\n")
    file(GENERATE
      OUTPUT "${CMAKE_BINARY_DIR}/${_python_module_path}/${_python_module_name}/${module}.py"
      CONTENT "${module_py}")
    install(
      FILES       "${CMAKE_BINARY_DIR}/${_python_module_path}/${_python_module_name}/${module}.py"
      DESTINATION "${CMAKE_INSTALL_PREFIX}/${_python_module_path}/${_python_module_name}")
  endif ()
else ()
  vtk_module_add_module(${module}
    SOURCES ${sources}
    HEADERS ${headers}
  )
endif ()

target_include_directories(${module} PRIVATE "${CMAKE_CURRENT_LIST_DIR}/../include")
target_link_libraries(${module} PRIVATE Eigen3::Eigen Boost::boost)

if (VTK_WRAP_PYTHON)
  # generate __init__.py that imports our module
  set(module_py "from __future__ import absolute_import\nfrom .${module} import *\n")
  file(GENERATE
    OUTPUT "${CMAKE_BINARY_DIR}/${_python_module_path}/${_python_module_name}/__init__.py"
    CONTENT "${module_py}")
  install(
    FILES       "${CMAKE_BINARY_DIR}/${_python_module_path}/${_python_module_name}/__init__.py"
    DESTINATION "${CMAKE_INSTALL_PREFIX}/${_python_module_path}/${_python_module_name}")
endif ()
