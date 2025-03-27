#!/bin/bash

# Create an array of the number of particles
particles=(1000 10000 100000 250000 500000)

# Number of iterations
iterations=100

# Create the logs directory if it doesn't exist
mkdir -p logs
mkdir -p logs/omp

#build folder
build_dir="../scripts/build_omp"

threads=(1 2 4 8 16 32 64)

# Run the simulation for each number of particles with varying number of threads
for particle in ${particles[@]}; do
    for thread in ${threads[@]}; do
        # Run the simulation 10 times
        for i in {1..10}; do
            # Run the simulation with the current number of particles and threads
            output=$(../scripts/run_omp.sh $thread $particle $iterations $build_dir)

            # Extract execution times
            seconds=$(echo "$output" | grep "Time: (seconds)" | awk '{print $NF}')
            milliseconds=$(echo "$output" | grep "Time: (milliseconds)" | awk '{print $NF}')
            microseconds=$(echo "$output" | grep "Time: (microseconds)" | awk '{print $NF}')

            # Ensure output directory exists
            dir="out/out_omp"
            mkdir -p "$dir"

            # Save the extracted times to a separate file
            echo "$seconds, $milliseconds, $microseconds" > "${dir}/times_${particle}_${thread}_${i}.txt"

            # Save the full output log in the logs directory
            echo "$output" > "logs/omp/output_${particle}_${thread}_${i}.txt"
        done
    done
done
