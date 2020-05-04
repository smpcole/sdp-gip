#!/bin/bash

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module load matlab/2020a
matlab -nodisplay -nosplash -nodesktop -r "feasibility(A, B, Inf, true, 'std'); exit"
