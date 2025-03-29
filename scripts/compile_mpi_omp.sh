#!/bin/bash
#parse arguments for sdl
USE_SDL=OFF
if [ "$1" == "sdl" ]; then
    USE_SDL=ON
fi
cmake -D CMAKE_C_COMPILER=mpicc -D CMAKE_CXX_COMPILER=mpicxx -DCMAKE_BUILD_TYPE=Release -S .. -B build_mpi_omp -DUSE_SDL=$USE_SDL -DUSE_MPI=ON -DUSE_OPENMP=ON
cmake --build ./build_mpi_omp
