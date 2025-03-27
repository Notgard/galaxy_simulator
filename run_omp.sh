#!/bin/bash

# Check if the correct number of arguments is provided
if [[ "$#" -lt 3 || "$#" -gt 4 ]]; then
    echo "Usage: $0 <number_of_threads> <number_of_particles> <number_of_iterations>"
    exit 1
fi

exec=./build/particleSimulation

# Assign arguments to variables
number_of_threads=$1
number_of_particles=$2
number_of_iterations=$3

OMP_NUM_THREADS=$number_of_threads $exec $number_of_particles $number_of_iterations