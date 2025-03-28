#!/usr/bin/env bash
#SBATCH --account="r250059"
#SBATCH --time=24:00:00
#SBATCH --mem=40G
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --job-name "bench_omp"
#SBATCH --comment "100000 100"
#SBATCH --error=../out/omp_job/job.err
#SBATCH --output=../out/omp_job/job.out

source ../load_env.sh

./bench_omp.sh
