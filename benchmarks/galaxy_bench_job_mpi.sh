#!/usr/bin/env bash
#SBATCH --account="r250059"
#SBATCH --time=12:00:00
#SBATCH --mem=30G
#SBATCH --constraint=x64cpu
#SBATCH --ntasks-per-node=32 #MPI
#SBATCH --ntasks=64 #MPI
#SBATCH --cpus-per-task=4 #OpenMP
#SBATCH --job-name "bench_mpi"
#SBATCH --comment "100000 100"
#SBATCH --error=out/mpi_job/job.err
#SBATCH --output=out/mpi_job/job.out

source ../load_env.sh

./bench_mpi.sh
