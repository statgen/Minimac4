#!/usr/bin/env cmake -P

set(ARGS)
foreach(i RANGE 4 ${CMAKE_ARGC})
    list(APPEND ARGS ${CMAKE_ARGV${i}})
endforeach()

set(_PREFIX ${CMAKE_ARGV3})

# Make sure this is in the module path
list(APPEND CMAKE_MODULE_PATH ${CMAKEGET_MODULE_PATH})
include(${CMAKE_CURRENT_LIST_DIR}/CMakeGet.cmake)

get_filename_component(PREFIX ${_PREFIX} ABSOLUTE)
# Install recipes
#cmake_get(pfultz2/cget-recipes PREFIX ${PREFIX} CMAKE_ARGS ${ARGS})
cmake_get_from(${CMAKE_CURRENT_LIST_DIR}/requirements.txt PREFIX ${PREFIX} CMAKE_ARGS ${ARGS})
