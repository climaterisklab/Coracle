#!/bin/bash
#SBATCH --job-name=reef_counts_v2
#SBATCH --account=bsc32
#SBATCH --qos=bsc_es
#SBATCH --array=0-11
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/scratch/bsc32/bsc963510/coracle_fig3/v4_area/logs/%x_%A_%a.out
# 84 seeds over 12 array tasks (7 seeds each). Afterwards: python reef_counts_v2.py --combine
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
PY=/gpfs/scratch/bsc32/bsc963510/conda/envs/coracle_plot/bin/python
SEEDS=$(seq $((SLURM_ARRAY_TASK_ID * 7)) $((SLURM_ARRAY_TASK_ID * 7 + 6)))
$PY /gpfs/scratch/bsc32/bsc963510/coracle_fig3/v4_area/reef_counts_v2.py $SEEDS
