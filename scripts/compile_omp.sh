#!/bin/bash
#parse arguments for sdl
USE_SDL=OFF
if [ "$1" == "sdl" ]; then
    USE_SDL=ON
fi
cmake -DCMAKE_BUILD_TYPE=Release -D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++ -S .. -B build_omp -DUSE_SDL=$USE_SDL -DUSE_MPI=OFF -DUSE_OPENMP=ON
cmake --build ./build_omp