#!/bin/bash

# Create an array of the number of particles
particles=(1000 10000 100000 250000 500000)

# Number of iterations
iterations=100

# Create the logs directory if it doesn't exist
mkdir -p logs
mkdir -p logs/seq

# Run the simulation for each number of particles
for particle in ${particles[@]}; do
    # Run the simulation 10 times
    for i in {1..10}; do
        # Run the simulation with the current number of particles and iterations
        output=$(../scripts/build/particleSimulation $particle 100)

        # Extract execution times
        seconds=$(echo "$output" | grep "Time: (seconds)" | awk '{print $NF}')
        milliseconds=$(echo "$output" | grep "Time: (milliseconds)" | awk '{print $NF}')
        microseconds=$(echo "$output" | grep "Time: (microseconds)" | awk '{print $NF}')

        # Save the extracted times to a separate file
        echo "$seconds, $milliseconds, $microseconds" > "out/out_seq/times_${particle}_${i}.txt"

        # Save the full output log in the logs directory
        echo "$output" > "logs/seq/output_${particle}_${i}.txt"
    done
done
