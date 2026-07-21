#!/bin/bash
#SBATCH --job-name=setup_linear_lat_conley_season
#SBATCH --output=/gpfs/scratch/bsc32/bsc963510/logs/setup_linear_lat_conley_season_%j.out
#SBATCH --error=/gpfs/scratch/bsc32/bsc963510/logs/setup_linear_lat_conley_season_%j.err
#SBATCH --cpus-per-task=20
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --time=12:00:00

echo "Setup linear_lat_conley_season on $(hostname) at $(date)"
module purge
module load CONDA-FORGE/miniforge3-23.3.1-1
conda activate coracle
ulimit -v unlimited
Rscript /gpfs/scratch/bsc32/bsc963510/setup_linear_lat_conley_season.R
