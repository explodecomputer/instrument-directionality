#!/bin/bash

#SBATCH --job-name=agg
#SBATCH --nodes=1 --cpus-per-task=1 --time=0-6:00:00
#SBATCH --partition=mrcieu
#SBATCH --output=job_reports/slurm-%A_%a.out
#SBATCH --mem=20G


set -e

echo "Running on ${HOSTNAME}"
date

cd ${HOME}/mr-eve/instrument-directionality/scripts

time Rscript sim-aggregate.r ../results/sim3

date
