#!/bin/bash
#parse arguments for sdl
USE_SDL=OFF
if [ "$1" == "sdl" ]; then
    USE_SDL=ON
fi
cmake -DCMAKE_BUILD_TYPE=Release -D -S . -B build -DUSE_SDL=$USE_SDL -DUSE_MPI=ON -DUSE_OPENMP=OFF
cmake --build build