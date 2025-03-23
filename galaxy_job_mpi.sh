#!/usr/bin/env bash
#SBATCH --account="r250059"
#SBATCH --time=1-02:30:00
#SBATCH --mem=30G
#SBATCH --constraint=x64cpu
#SBATCH --ntasks-per-node=32 #MPI
#SBATCH --ntasks=64 #MPI
#SBATCH --cpus-per-task=4 #OpenMP
#SBATCH --job-name "galaxy_mpi"
#SBATCH --comment "100000 100"
#SBATCH --error=out/mpi_job/job.err
#SBATCH --output=out/mpi_job/job.out

romeo_load_x64cpu_env
spack load openmpi@5.0.5/rvuhwou

./run_mpi.sh 64 32 100000 100
