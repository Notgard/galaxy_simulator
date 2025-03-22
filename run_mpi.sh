#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <number_of_process> <number_of_threads> <number_of_particles> <number_of_iterations>"
    exit 1
fi

# Assign arguments to variables
number_of_process=$1
number_of_threads=$2
number_of_particles=$3
number_of_iterations=$4

# Execute the command
mpiexec -n $number_of_process -x OMP_NUM_THREADS=$number_of_threads ./build/particleSimulation $number_of_particles $number_of_iterations