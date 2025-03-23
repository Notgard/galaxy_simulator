#!/usr/bin/env bash
#SBATCH --account="r250059"
#SBATCH --time=1-02:30:00
#SBATCH --mem=30G
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name "galaxy_seq"
#SBATCH --comment "100000 100"
#SBATCH --error=out/seq_job/job.err
#SBATCH --output=out/seq_job/job.out

romeo_load_x64cpu_env

./build/particleSimulation 100000 100
