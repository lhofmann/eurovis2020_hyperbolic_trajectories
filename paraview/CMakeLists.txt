set(plugin             HyperbolicTrajectoriesPlugin)
set(plugin_version     "1.0")
set(server_manager_xml ApproximateDHT.xml)
set(modules            HyperbolicTrajectories)

if (_legacy_build_system)
  set(_headers)
  foreach (module ${modules})
    get_target_property(_source_dir ${module} SOURCE_DIR)
    get_target_property(_headers_relative ${module} HEADERS)
    foreach (header ${_headers_relative})
      if (IS_ABSOLUTE ${header})
        set(_headers ${_headers} "${header}")
      else ()
        set(_headers ${_headers} "${_source_dir}/${header}")
      endif ()
    endforeach ()
  endforeach ()

  add_paraview_plugin(${plugin} ${plugin_version}
      SERVER_MANAGER_XML ${server_manager_xml}
      SERVER_MANAGER_SOURCES ${_headers})
  target_link_libraries(${plugin} PRIVATE ${modules})
  set_target_properties(${plugin} PROPERTIES PREFIX "")  
  set_target_properties(${plugin} PROPERTIES INSTALL_RPATH "${_paraview_plugin_rpath}")
  set_target_properties(${plugin} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${_paraview_plugin_path})
  
  install(TARGETS ${plugin} DESTINATION ${CMAKE_INSTALL_PREFIX}/${_paraview_plugin_path})
else ()
  paraview_add_plugin(${plugin}
    VERSION ${plugin_version}
    MODULES ${modules}
    SERVER_MANAGER_XML ${server_manager_xml}
  )
  set_target_properties(${plugin} PROPERTIES INSTALL_RPATH "${_paraview_plugin_rpath}")
  set_target_properties(${plugin} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${_paraview_plugin_path})

  # paraview_add_plugin installs the plugin into bin/plugin/plugin.so
  # in order for ParaView to find the plugin, we create a symlink bin/plugin.so to bin/plugin/plugin.so
  install(CODE "execute_process( \
      COMMAND ${CMAKE_COMMAND} -E create_symlink \
      ${plugin}/${plugin}.so \
      ${CMAKE_INSTALL_PREFIX}/${_paraview_plugin_path}/${plugin}.so \
      )"
  )
endif ()
