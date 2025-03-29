#!/usr/bin/env bash
#SBATCH --account="r250059"
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name "bench_seq"
#SBATCH --error=../out/seq_job/job.err
#SBATCH --output=../out/seq_job/job.out

source ../load_env.sh

./bench_seq.sh
