#!/usr/bin/env bash
#SBATCH --account="r250059"
#SBATCH --time=24:00:00
#SBATCH --mem=40G
#SBATCH --constraint=x64cpu
#SBATCH --ntasks-per-node=64 #MPI
#SBATCH --ntasks=128 #MPI
#SBATCH --cpus-per-task=2 #OpenMP
#SBATCH --job-name "bench_mpi_omp_128"
#SBATCH --comment "100000 100"
#SBATCH --error=../out/mpi_omp_job_128/job.err
#SBATCH --output=../out/mpi_omp_job_128/job.out

source ../load_env.sh

#./bench_mpi_omp.sh
./bench_mpi_omp_128.sh
