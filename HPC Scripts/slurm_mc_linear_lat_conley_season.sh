#!/bin/bash
#SBATCH --job-name=mc_linear_lat_conley_season
#SBATCH --output=/gpfs/scratch/bsc32/bsc963510/logs/mc_linear_lat_conley_season_%A_%a.out
#SBATCH --error=/gpfs/scratch/bsc32/bsc963510/logs/mc_linear_lat_conley_season_%A_%a.err
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --array=1-84

echo "MC batch $SLURM_ARRAY_TASK_ID for linear_lat_conley_season on $(hostname) at $(date)"
module purge
module load CONDA-FORGE/miniforge3-23.3.1-1
conda activate coracle
ulimit -v unlimited
Rscript /gpfs/scratch/bsc32/bsc963510/mc_batch_linear_lat_conley_season.R
