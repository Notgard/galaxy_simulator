#!/bin/bash

# Create an array of the number of particles
particles=(1000 10000 100000 250000 500000)

# Number of iterations
iterations=100

# Number of processes
processes=(4 8 16 32 64 128)

# Number of threads
threads=1

#build folder
build_dir="../scripts/build_mpi"

# MPI versions
mpi_versions=(1 2)

mkdir -p logs
mkdir -p logs/mpi

# Run the simulation for each number of particles
for mpi_version in ${mpi_versions[@]}; do
    for particle in ${particles[@]}; do
        # Run the simulation for each number of processes
        for process in ${processes[@]}; do
            # Run the simulation 10 times
            for i in {1..10}; do
                # Run the simulation with the current number of particles and processes
                output=$(../scripts/run_mpi.sh $process $threads $particle $iterations $mpi_version $build_dir)

                # Extract execution times
                seconds=$(echo "$output" | grep "Time: (seconds)" | awk '{print $NF}')
                milliseconds=$(echo "$output" | grep "Time: (milliseconds)" | awk '{print $NF}')
                microseconds=$(echo "$output" | grep "Time: (microseconds)" | awk '{print $NF}')

                # Ensure output directory exists
                if [ $mpi_version -eq 1 ]; then
                    dir="out/out_mpi/threads_1/sub_tree"
                else
                    dir="out/out_mpi/threads_1/chunk"
                fi
                mkdir -p "$dir"

                # Save the extracted times to a separate file
                echo "$seconds, $milliseconds, $microseconds" > "${dir}/times_${particle}_${process}_${i}.txt"

                # Save the full output log in the logs directory
                echo "$output" > "logs/mpi/output_${particle}_${process}_${i}.txt"
            done
        done
    done
done
