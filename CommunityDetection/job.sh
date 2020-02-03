#!/bin/bash

#SBATCH --time=3:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=8192M   # memory per CPU core
#SBATCH --mail-user=michael.bradshawiii@colorado.edu   # email address
#SBATCH --mail-type=FAIL

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
Rscript beckett_community_detection.R $1"a"
Rscript beckett_community_detection.R $1"b"
Rscript beckett_community_detection.R $1"c"
Rscript beckett_community_detection.R $1"d"
Rscript beckett_community_detection.R $1"e"
Rscript beckett_community_detection.R $1"f"
Rscript beckett_community_detection.R $1"g"
Rscript beckett_community_detection.R $1"h"
Rscript beckett_community_detection.R $1"i"
Rscript beckett_community_detection.R $1"j"








