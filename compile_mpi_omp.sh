#!/bin/bash
#parse arguments for sdl
USE_SDL=OFF
if [ "$1" == "sdl" ]; then
    USE_SDL=ON
fi
cmake -DCMAKE_EXE_LINKER_FLAGS=-pg -DCMAKE_SHARED_LINKER_FLAGS=-pg -S . -B build -DUSE_SDL=$USE_SDL -DUSE_MPI=ON -DUSE_OPENMP=ON
cmake --build build