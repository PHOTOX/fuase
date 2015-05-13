#!/bin/bash

# ------------------------------------------------------------
# This is a sample script to show how to submit TurboMole 6.2
# job to SGE queueing system for parallel run.
# ------------------------------------------------------------

# After modyfying the program call procedure submit the script
# using the following command:
#
# qsub -q <queue> -pe hpmpi <n> <scriptname>
#
# Replace <queue> with the queue name,  <n> with the number of
# processors, and <scriptname> with the name of this script.

# ------------------------------------------------------------
# DO NOT MODIFY THE FOLLOWING LINES!! THE NUMBER OF PROCESSORS
# WILL BE SET ON THE COMMAND LINE DURING SUBMISSION.
# ------------------------------------------------------------

#$ -cwd
#$ -e .
#$ -o .
#$ -pe mpi 1
# $ -M <name>@<provider>.com
#$ -m eas
#$ -q aq

# . /home/uhlig/build/Python/2.7.9/env.sh
# . /home/uhlig/build/numpy/1.9.1/env.sh

. /home/uhlig/prog/fuase.trunk/ase/env.sh

. /home/uhlig/build/openmpi/1.6.5/env.sh
export PATH=/home/uhlig/build/orca/3.0.2:$PATH

# ------------------------------------------------------------

touch jobinfo_$(date +"%F_%T")_${HOSTNAME}_${QUEUE}_${JOB_ID}

JOBDIR=${PWD}
SCRDIR=/scratch/${USER}/TM_${QUEUE}_${JOB_ID}

mkdir -p ${SCRDIR}
cp -a . ${SCRDIR}
cd ${SCRDIR}

# ------------------------------------------------------------
# MODIFY HERE - CALL PROGRAMS, REPLACE "output" WITH
# THE DESIRED NAME OF OUTPUT FILE
# ------------------------------------------------------------

./run-ase-script.py

# ------------------------------------------------------------
# DO NOT MODIFY THE FOLLOWING LINES!!
# ------------------------------------------------------------

cd ${JOBDIR}

if cp -a ${SCRDIR}/* . ;  then
  rm -rf ${SCRDIR}
else
  touch cp_fail_$(date +"%F_%T")_${HOSTNAME}_${QUEUE}_${JOB_ID}
fi

# ------------------------------------------------------------

