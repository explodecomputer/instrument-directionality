#!/bin/bash

#SBATCH --job-name=sim4
#SBATCH --array=1-500
#SBATCH --nodes=1 --cpus-per-task=1 --time=0-12:00:00
#SBATCH --partition=mrcieu
#SBATCH --output=job_reports/slurm-%A_%a.out
#SBATCH --mem=24G

echo "Running on ${HOSTNAME}"
module add languages/r/3.4.4

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}

sims=100

cd ${HOME}/mr-eve/instrument-directionality/scripts
mkdir -p ../results/sim4
Rscript \
        sim-run2.r \
        ${i} \
        ${sims} \
        ../results/sim4/simulate4_${i}-metrics.rdata

