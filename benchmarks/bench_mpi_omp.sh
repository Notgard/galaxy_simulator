#!/bin/bash

# create an array of the number of particles
particles=(1000 10000 100000 250000 500000)

#first benchmark with 4 threads
#number of iterations
iterations=100

#number of MPI processes
processes=(4 8 16 32 64)

#number of OMP threads
threads=4

#MPI version
mpi_versions=(1 2)

#run the simulation for each number of particles
for mpi_version in ${mpi_versions[@]}
do
    for particle in ${particles[@]}
    do
        #run the simulation for each number of processes
        for process in ${processes[@]}
        do
            #run the simulation 10 times
            for i in {1..10}
            do
                #run the simulation with the current number of particles and processes
                ../run_mpi.sh $process $threads $particle $iterations $mpi_version  > output.txt

                
                #put the output file in a folder
                if(mpi_version == 1)
                then
                    mv output.txt out/out_mpi/threads_4/sub_tree/output_${particle}_${process}_${i}.txt
                else
                    mv output.txt out/out_mpi/threads_4/chunk/output_${particle}_${process}_${i}.txt

            done
        done
    done
done
#second benchmark with 2 threads

#number of processes
processes=(4 8 16 32 64 128)

#number of threads
threads=2

#run the simulation for each number of particles
for mpi_version in ${mpi_versions[@]}
do
    for particle in ${particles[@]}
    do
        #run the simulation for each number of processes
        for process in ${processes[@]}
        do
            #run the simulation 10 times
            for i in {1..10}
            do
                #run the simulation with the current number of particles and processes
                ../run_mpi.sh $process $threads $particle $iterations > output.txt

                 #put the output file in a folder
                if(mpi_version == 1)
                then
                    mv output.txt out/out_mpi/threads_2/sub_tree/output_${particle}_${process}_${i}.txt
                else
                    mv output.txt out/out_mpi/threads_2/chunk/output_${particle}_${process}_${i}.txt

            done
        done
    done
done