#!/bin/bash

#SBATCH --job-name=sims
#SBATCH --array=501-1000
#SBATCH --nodes=1 --cpus-per-task=1 --time=0-12:00:00
#SBATCH --partition=mrcieu
#SBATCH --output=job_reports/slurm-%A_%a.out
#SBATCH --mem=12G


set -e

echo "Running on ${HOSTNAME}"
date

if [ -n "${1}" ]; then
	echo "${1}"
	SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}

cd ${HOME}/mr-eve/instrument-directionality/scripts

time Rscript \
	sim-metrics.r \
	../results/sim3/simulate3_${i}.rdata \
	../results/sim3/simulate3_${i}-metrics.rdata \
	${i}

date
