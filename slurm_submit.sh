#!/bin/bash
#SBATCH --partition general
#SBATCH --array 1-20
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 0-12:00
#SBATCH --output=./output_logs/output-%A_%a.out

# Load conda environment
module load miniconda3
source activate /home/ayn6k/.conda/envs/spindle

echo "Running: ${SLURM_JOB_NAME} (Task $SLURM_ARRAY_TASK_ID)"
echo "Arguments: $@"

# Get today's date and setup directories  
today_date=$(date +'%Y-%m-%d')
target_directory=$(pwd)
date_folder="$target_directory/output/$today_date"
ensemble_folder="$date_folder/${SLURM_JOB_NAME}"

# Create directory structure
mkdir -p "$date_folder"
mkdir -p "$ensemble_folder"
mkdir -p "$ensemble_folder/${SLURM_JOB_NAME}_images"
mkdir -p "$ensemble_folder/${SLURM_JOB_NAME}_stats"

# Pass ALL arguments directly to Python (the magic!)
python3 spindle.py "$@"

echo "Simulation complete for task $SLURM_ARRAY_TASK_ID"
