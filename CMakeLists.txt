cmake_minimum_required(VERSION 3.2)
project(metis_partition)
set(CMAKE_CXX_STANDARD 20)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

#get kahip deploy path
if (NOT DEFINED KAHIP_PATH)
    set(KAHIP_PATH $ENV{HOME}/local)
endif()
include_directories("${KAHIP_PATH}")

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
        "${KAHIP_PATH}/libparhip.a"
        "${KAHIP_PATH}/parallel/libkahip.a"
)

add_executable(convert src/convert.cpp src/io.cpp src/io.h src/quicksort.h src/quicksort.cpp src/typedef.h src/typedef.cpp)
target_link_libraries(convert PRIVATE
        ${MPI_C_LIBRARIES}
        OpenMP::OpenMP_CXX
        "${KAHIP_PATH}/libparhip.a"
        "${KAHIP_PATH}/parallel/libkahip.a"
)

target_include_directories(metis_partition PRIVATE ${MPI_C_INCLUDE_PATH})
