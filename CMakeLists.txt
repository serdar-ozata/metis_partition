cmake_minimum_required(VERSION 3.2)
project(metis_partition)
set(CMAKE_CXX_STANDARD 20)

find_package(MPI REQUIRED)

# Default library paths
if (NOT DEFINED PARMETIS_PATH)
    set(PARMETIS_PATH $ENV{HOME}/local)
endif()

if (NOT DEFINED METIS_PATH)
    set(METIS_PATH $ENV{HOME}/local)
endif()

if (NOT DEFINED GKLIB_PATH)
    set(GKLIB_PATH $ENV{HOME}/local)
endif()

include_directories(${PARMETIS_PATH}/include)
include_directories(${METIS_PATH}/include)
include_directories(${GKLIB_PATH}/include)

# Locate libraries
find_library(PARMETIS_SHARED NAMES parmetis PATHS "${PARMETIS_PATH}/lib" NO_DEFAULT_PATH)
find_library(PARMETIS_STATIC NAMES libparmetis.a PATHS "${PARMETIS_PATH}/lib" NO_DEFAULT_PATH)

find_library(METIS_SHARED NAMES metis PATHS "${METIS_PATH}/lib" NO_DEFAULT_PATH)
find_library(METIS_STATIC NAMES libmetis.a PATHS "${METIS_PATH}/lib" NO_DEFAULT_PATH)

find_library(GKLIB_SHARED NAMES GKlib PATHS "${GKLIB_PATH}/lib" NO_DEFAULT_PATH)
find_library(GKLIB_STATIC NAMES libGKlib.a PATHS "${GKLIB_PATH}/lib" NO_DEFAULT_PATH)

if (PARMETIS_SHARED)
    set(PARMETIS_LIB "${PARMETIS_SHARED}")
elseif (PARMETIS_STATIC)
    set(PARMETIS_LIB "${PARMETIS_STATIC}")
else()
    message(FATAL_ERROR "ParMETIS library not found in ${PARMETIS_PATH}/lib")
endif()

if (METIS_SHARED)
    set(METIS_LIB "${METIS_SHARED}")
elseif (METIS_STATIC)
    set(METIS_LIB "${METIS_STATIC}")
else()
    message(FATAL_ERROR "METIS library not found in ${METIS_PATH}/lib")
endif()

if (GKLIB_SHARED)
    set(GKLIB_LIB "${GKLIB_SHARED}")
elseif (GKLIB_STATIC)
    set(GKLIB_LIB "${GKLIB_STATIC}")
else()
    message(FATAL_ERROR "GKlib library not found in ${GKLIB_PATH}/lib")
endif()

add_executable(metis_partition src/main.cpp
        src/io.cpp
        src/io.h
        src/typedef.h
        src/partition.cpp
        src/partition.h
        src/typedef.cpp
        src/quicksort.cpp
        src/quicksort.h)

target_link_libraries(metis_partition PRIVATE
        ${MPI_C_LIBRARIES}
        "${PARMETIS_LIB}"
        "${METIS_LIB}"
        "${GKLIB_LIB}"
)

target_include_directories(metis_partition PRIVATE ${MPI_C_INCLUDE_PATH})
