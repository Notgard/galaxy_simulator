#!/usr/bin/env bash
romeo_load_x64cpu_env
spack load cmake
spack load aocc@5.0.0
spack load openmpi@5.0.5 %aocc
