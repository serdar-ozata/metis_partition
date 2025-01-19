cmake_minimum_required(VERSION 3.2)
project(metis_partition)
set(CMAKE_CXX_STANDARD 20)

find_package(MPI REQUIRED)

if (NOT DEFINED PARMETIS_PATH)
    set(PARMETIS_PATH $ENV{HOME}/local)
endif()
include_directories(${PARMETIS_PATH}/include)

if (NOT DEFINED METIS_PATH)
    set(METIS_PATH $ENV{HOME}/local)
endif()

if (NOT DEFINED GKLIB_PATH)
    set(GKLIB_PATH $ENV{HOME}/local)
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

#target_link_libraries(metis_partition PRIVATE ${MPI_C_LIBRARIES} $ENV{HOME}/local/lib/libparmetis.a /lib/libmetis.so /lib/libGKlib.so)
target_link_libraries(metis_partition PRIVATE ${MPI_C_LIBRARIES} ${PARMETIS_PATH}/lib/libparmetis.a ${METIS_PATH}/lib/libmetis.so ${GKLIB_PATH}/lib/libGKlib.so)
target_include_directories(metis_partition PRIVATE ${MPI_C_INCLUDE_PATH})