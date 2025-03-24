#!/bin/bash

# Check if the correct number of arguments is provided
if [[ "$#" -lt 4 || "$#" -gt 5 ]]; then
    echo "Usage: $0 <number_of_process> <number_of_threads> <number_of_particles> <number_of_iterations> [sdl]"
    exit 1
fi

use_sdl=""
if [ "$5" == "sdl" ]; then
    use_sdl="sdl"
fi

exec=./build/particleSimulationMPI

# Assign arguments to variables
number_of_process=$1
number_of_threads=$2
number_of_particles=$3
number_of_iterations=$4

# Execute the command
mpiexec -n $number_of_process -x OMP_NUM_THREADS=$number_of_threads $exec $number_of_particles $number_of_iterations $use_sdl