#!/usr/bin/env bash
#SBATCH --account="r250059"
#SBATCH --time=10:00:00
#SBATCH --mem=30G
#SBATCH --constraint=x64cpu
#SBATCH --ntasks-per-node=32 #MPI
#SBATCH --ntasks=64 #MPI
#SBATCH --cpus-per-task=4 #OpenMP
#SBATCH --job-name "galaxy_mpi_omp"
#SBATCH --comment "100000 100"
#SBATCH --error=out/mpi_omp_job/job.err
#SBATCH --output=out/mpi_omp_job/job.out

romeo_load_x64cpu_env
spack load openmpi@5.0.5/rvuhwou

./run_mpi.sh 64 4 100000 100
