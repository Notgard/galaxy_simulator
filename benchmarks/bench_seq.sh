#!/bin/bash

# create an array of the number of particles
particles=(1000 10000 100000 250000 500000)

#number of iterations
iterations=100

#run the simulation for each number of particles
for particle in ${particles[@]}
do
    #run the simulation 10 times
    for i in {1..10}
    do
        #run the simulation with the current number of particles and iterations
        ../build/particleSimulation $particle 100 > output.txt

        #put the output file in a folder 
        mv output.txt out/out_seq/output_${particle}_${i}.txt

    done
done