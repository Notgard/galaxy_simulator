#!/bin/bash

# Check if the correct number of arguments is provided
if [[ "$#" -lt 4 || "$#" -gt 5 ]]; then
    echo "Usage: $0 <number_of_threads> <number_of_particles> <number_of_iterations> <build_dir>"
    exit 1
fi

# Assign arguments to variables
number_of_threads=$1
number_of_particles=$2
number_of_iterations=$3
build_dir=$4

exec=$build_dir/particleSimulation

OMP_PLACES=cores OMP_PROC_BIND=close OMP_NUM_THREADS=$number_of_threads $exec $number_of_particles $number_of_iterations