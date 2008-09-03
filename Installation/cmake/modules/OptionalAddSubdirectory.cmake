macro( add_subdirectory_if cond dir)

  if ( ${cond} )
    message( STATUS "Configuring ${dir}. Set ${cond} to FALSE to unselect it." )
    add_subdirectory( ${dir} ${ARGN} )
  endif()

endmacro()

macro( optional_add_subdirectory dir def)
  set( WITH_${dir}_ENV $ENV{WITH_${dir}} )
  if ( NOT ${WITH_${dir}_ENV} STREQUAL "" )
    message ( STATUS "WITH_${dir}_ENV given as enviroment variable: ${WITH_${dir}_ENV}" )
    set( WITH_${dir} ${WITH_${dir}_ENV} CACHE BOOL "Select ${dir} package." FORCE )
  endif()
  option( WITH_${dir} "Select ${dir} package." ${def} )
  add_subdirectory_if( WITH_${dir} ${dir} ${ARGN} )
endmacro()

