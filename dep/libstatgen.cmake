cmake_minimum_required(VERSION 3.2)
project(libStatGen VERSION 1.0.0)

#execute_process(COMMAND ./configure --prefix=${CMAKE_INSTALL_PREFIX} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_custom_target(libStatGen ALL COMMAND make WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMENT "Builing libStatGen ...")

file(GLOB_RECURSE LSG_HEADER_LIST "bam/*.h" "fastq/*.h" "general/*.h" "glf/*.h" "samtools/*.h" "vcf/*.h")
message("HEADERS: ${LSG_HEADER_LIST}")
install(FILES ${LSG_HEADER_LIST} DESTINATION include)

if (BUILD_SHARED_LIBS)
    install(FILES ${CMAKE_SHARED_LIBRARY_PREFIX}StatGen${CMAKE_SHARED_LIBRARY_SUFFIX} DESTINATION lib)
else()
    install(FILES ${CMAKE_STATIC_LIBRARY_PREFIX}StatGen${CMAKE_STATIC_LIBRARY_SUFFIX} DESTINATION lib)
endif()