#!/bin/bash

#PBS -N case_collider
#PBS -o job_reports/casecollider-output
#PBS -e job_reports/casecollider-error
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

Rscript \
	power_calc.r \
	${i} \
	${splits} \
	../results/power_calc_${i}.rdata
