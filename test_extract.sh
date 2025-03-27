#!/bin/bash

# Run the C++ program and capture output
output=$(./build/particleSimulation 100 1)

# Extract execution times
seconds=$(echo "$output" | grep "Time: (seconds)" | awk '{print $NF}')
milliseconds=$(echo "$output" | grep "Time: (milliseconds)" | awk '{print $NF}')
microseconds=$(echo "$output" | grep "Time: (microseconds)" | awk '{print $NF}')

# Print extracted values
echo "Seconds: $seconds"
echo "Milliseconds: $milliseconds"
echo "Microseconds: $microseconds"

