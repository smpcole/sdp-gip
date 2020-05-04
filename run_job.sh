#!/bin/bash

if [ ! -z $SLURM_CPUS_PER_TASK ]; then
    export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
fi

NUM_NONNEG=Inf
PSD=true
BASIS=std

if [ $# -ge 1 ]; then
    NUM_NONNEG="$1"
fi
if [ $# -ge 2 ]; then
    PSD="$2"
fi
if [ $# -eq 3 ]; then
    BASIS="$3"
fi
if [ $# -gt 3 ]; then
    echo Usage: ./run_job.sh [num_nonneg] [psd] [basis]
    exit 1
fi

CMD="feasibility(A, B, $NUM_NONNEG, $PSD, '$BASIS'); exit"

module load matlab/2020a
matlab -nodisplay -nosplash -nodesktop -r "$CMD"
