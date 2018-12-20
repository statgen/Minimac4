cmake_minimum_required(VERSION 3.2)
project(minimac4 VERSION 1.0.1)
include(CMakePackageConfigHelpers)

set(CMAKE_CXX_STANDARD 11)
execute_process(COMMAND date OUTPUT_VARIABLE DATE OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND whoami OUTPUT_VARIABLE USER OUTPUT_STRIP_TRAILING_WHITESPACE)

add_definitions(-DVERSION="${PROJECT_VERSION}" -DUSER="${USER}" -DDATE="${DATE}")

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;${CMAKE_FIND_LIBRARY_SUFFIXES}")

find_package(savvy REQUIRED)
find_library(STATGEN_LIBRARY StatGen)

add_executable(minimac4
        src/Analysis.cpp
        src/Analysis.h
        src/AnalysisChunks.cpp
        src/DosageData.cpp
        src/DosageData.h
        src/HaplotypeSet.cpp
        src/HaplotypeSet.h
        src/Imputation.cpp
        src/Imputation.h
        src/ImputationStatistics.cpp
        src/ImputationStatistics.h
        src/Main.cpp
        src/MarkovModel.cpp
        src/MarkovModel.h
        src/MarkovParameters.cpp
        src/MarkovParameters.h
        src/MyVariables.h
        src/Unique.cpp
        src/Unique.h
        src/Estimation.cpp
        src/Estimation.h)

target_link_libraries(minimac4 ${STATGEN_LIBRARY} savvy ${ZLIB_LIBRARIES})

add_custom_target(manuals
                  COMMAND help2man --output "${CMAKE_BINARY_DIR}/minimac4.1" "${CMAKE_BINARY_DIR}/minimac4")

install(TARGETS minimac4 RUNTIME DESTINATION bin)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/minimac4.1
        DESTINATION share/man/man1
        OPTIONAL)

set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
include(CPack)
