#!/usr/bin/env bash
#SBATCH --account="r250059"
#SBATCH --time=24:00:00
#SBATCH --mem=40G
#SBATCH --constraint=x64cpu
#SBATCH --ntasks-per-node=32 #MPI
#SBATCH --ntasks=64 #MPI
#SBATCH --cpus-per-task=4 #OpenMP
#SBATCH --job-name "bench_mpi_omp"
#SBATCH --comment "100000 100"
#SBATCH --error=../out/mpi_omp_job/job.err
#SBATCH --output=../out/mpi_omp_job/job.out

source ../load_env.sh

./bench_mpi_omp.sh
