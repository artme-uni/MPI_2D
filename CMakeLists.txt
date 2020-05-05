cmake_minimum_required(VERSION 3.15)
project(MPI_2D)

SET(CMAKE_CXX_COMPILER /usr//local/bin/mpicxx)

include_directories(/usr/local/opt/mpich/include)
link_directories(/usr/local/opt/mpich/lib)

add_executable(MPI_2D main.cpp mat_multiplication.cpp mat_multiplication.h)
target_link_libraries(MPI_2D mpi)

