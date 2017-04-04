cmake_minimum_required(VERSION 3.6)
project(open_mpi)

set(CMAKE_CXX_STANDARD 11)

include_directories(/usr/local/openmpi/include)

set(CMAKE_C_COMPILER /usr/local/open-mpi/bin/mpicc)
set(CMAKE_CXX_COMPILER /usr/local/open-mpi/bin/mpicxx)


set(SOURCE_FILES page_rank.cpp test.cpp page_rank_single.cpp)
add_executable(open_mpi ${SOURCE_FILES})