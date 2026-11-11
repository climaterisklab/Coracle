#!/bin/bash
#SBATCH --job-name=mc_v2_mixed
#SBATCH --account=bsc32
#SBATCH --qos=bsc_es
#SBATCH --output=/gpfs/scratch/bsc32/bsc963510/logs/mc_v2_mixed_%A_%a.out
#SBATCH --error=/gpfs/scratch/bsc32/bsc963510/logs/mc_v2_mixed_%A_%a.err
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=1-84

echo "MC batch $SLURM_ARRAY_TASK_ID for v2_mixed on $(hostname) at $(date)"
module purge
module load CONDA-FORGE/miniforge3-23.3.1-1
conda activate coracle
ulimit -v unlimited
Rscript /gpfs/scratch/bsc32/bsc963510/mc_batch_v2_mixed.R
