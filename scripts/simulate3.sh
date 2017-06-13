#!/bin/bash

#PBS -N inst-dir
#PBS -o job_reports/inst-dir-output
#PBS -e job_reports/inst-dir-error
#PBS -l walltime=12:00:00
#PBS -t 1-100
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
	simulate3.r \
	${i} \
	${splits} \
	../results/simulate3_${i}.rdata
