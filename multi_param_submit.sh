#!/bin/bash

# Multi-parameter sweep - edit these arrays for any parameters!
n_astral_mts=(100)
pull_force=(5)
push_force=(0.0)
spindle_half_length=(0.45)

# Base parameters
CELL_TYPE="follicle_epithelial"

echo "Multi-parameter sweep:"
echo "  n_astral_mts: ${n_astral_mts[@]}"
echo "  pull_force: ${pull_force[@]}"
echo "  spindle_half_length: ${spindle_half_length[@]}"
echo ""

# Loop over parameter combinations
for n_mts in "${n_astral_mts[@]}"; do
  for pull in "${pull_force[@]}"; do  
    for push in "${push_force[@]}"; do
      for sl in "${spindle_half_length[@]}"; do
        # Create job name
        job_name="${CELL_TYPE}_MT_${n_mts}_pull_${pull}_push_${push}_SL_${sl}_test_1"
        
        # Build argument string (any parameters can go here!)
        args="--cell_type $CELL_TYPE --n_astral_mts $n_mts --pull_force $pull  --push_force $push"
        
        echo "Submitting: $job_name"
        
        # Submit job with all arguments passed through
        sbatch --job-name="$job_name" slurm_submit.sh $args
      done
    done
  done
done

echo ""
echo "All jobs submitted!"
echo "Total combinations: $((${#n_astral_mts[@]} * ${#pull_force[@]} * ${#spindle_half_length[@]}))"
