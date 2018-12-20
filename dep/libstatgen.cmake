cmake_minimum_required(VERSION 3.2)
project(libStatGen VERSION 1.0.0)

#execute_process(COMMAND ./configure --prefix=${CMAKE_INSTALL_PREFIX} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I${dir}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${dir}")
endforeach()
message("dirs: ${dirs}")
message("cxxflags: ${CMAKE_CXX_FLAGS}")

set(ENV{CFLAGS}  "${CMAKE_C_FLAGS}")
set(ENV{CXXFLAGS} "${CMAKE_CXX_FLAGS}")
set(ENV{LDFLAGS} "-L${CMAKE_PREFIX_PATH}/lib")

add_custom_target(libStatGen ALL COMMAND make CFLAGS="\"${CMAKE_C_FLAGS}\"" WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMENT "Builing libStatGen ...")

file(GLOB_RECURSE LSG_HEADER_LIST "bam/*.h" "fastq/*.h" "general/*.h" "glf/*.h" "samtools/*.h" "vcf/*.h")
install(FILES ${LSG_HEADER_LIST} DESTINATION include)

if (BUILD_SHARED_LIBS)
    install(FILES ${CMAKE_SHARED_LIBRARY_PREFIX}StatGen${CMAKE_SHARED_LIBRARY_SUFFIX} DESTINATION lib)
else()
    install(FILES ${CMAKE_STATIC_LIBRARY_PREFIX}StatGen${CMAKE_STATIC_LIBRARY_SUFFIX} DESTINATION lib)
endif()
