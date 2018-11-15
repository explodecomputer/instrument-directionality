#!/bin/bash

#PBS -N inst-dir
#PBS -o job_reports/sim3m-output
#PBS -e job_reports/sim3m-error
#PBS -l walltime=12:00:00
#PBS -t 401-500
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash

set -e

echo "Running on ${HOSTNAME}"

if [ -n "${1}" ]; then
	echo "${1}"
	PBS_ARRAYID=${1}
fi

i=${PBS_ARRAYID}
splits=100

cd ${HOME}/repo/instrument-directionality/scripts

Rscript \
	sim-metrics.r \
	../results/scratch/sim_${i}.rdata \
	../results/scratch/sim_${i}-metrics.rdata \
	${i}
