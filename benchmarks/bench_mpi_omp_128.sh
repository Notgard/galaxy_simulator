#!/bin/bash

# Create an array of the number of particles
particles=(1000 10000 100000 250000 500000)

# Number of iterations
iterations=100

#build folder
build_dir="../scripts/build_mpi_omp"

# MPI version
mpi_versions=(1 2)

# Number of processes
processes=(4 8 16 32 64 128)

# Number of threads
threads=2

# Create the logs directory if it doesn't exist
mkdir -p logs
mkdir -p logs/mpi_omp

# Run the simulation for each number of particles
for mpi_version in ${mpi_versions[@]}; do
    for particle in ${particles[@]}; do
        # Run the simulation for each number of processes
        for process in ${processes[@]}; do
            # Run the simulation 10 times
            for i in {1..10}; do
                # Run the simulation with the current number of particles and processes
                output=$(../run_mpi.sh $process $threads $particle $iterations $mpi_version)

                # Extract execution times
                seconds=$(echo "$output" | grep "Time: (seconds)" | awk '{print $NF}')
                milliseconds=$(echo "$output" | grep "Time: (milliseconds)" | awk '{print $NF}')
                microseconds=$(echo "$output" | grep "Time: (microseconds)" | awk '{print $NF}')

                # Ensure output directory exists
                if [ $mpi_version -eq 1 ]; then
                    dir="out/out_mpi_omp/threads_2/sub_tree"
                else
                    dir="out/out_mpi_omp/threads_2/chunk"
                fi
                mkdir -p "$dir"

                # Save the extracted times to a separate file
                echo "$seconds, $milliseconds, $microseconds" > "${dir}/times_${particle}_${process}_${i}.txt"

                # Save the full output log in the logs directory
                echo "$output" > "logs/mpi_omp/output_${particle}_${process}_${i}.txt"
            done
        done
    done
done