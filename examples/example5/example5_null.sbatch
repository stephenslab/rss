#!/bin/bash

#SBATCH --job-name=ibd2015
#SBATCH --output=ibd2015_%A_%a.out
#SBATCH --error=ibd2015_%A_%a.err
#SBATCH --array=1-3
#SBATCH --partition=sandyb
#SBATCH --time=1-12:00:00
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=10
#SBATCH --mem-per-cpu=3200 
#SBATCH --mail-type=ALL,ARRAY_TASKS  
#SBATCH --mail-user=xiangzhu.nku@gmail.com

module load matlab/2014b

# Print this sub-job's task ID
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

# Copy the workhorse script to this sub-job
cp example5_null.m case_"$SLURM_ARRAY_TASK_ID".m

# Add two lines in the sub-job script
echo "case_id = $SLURM_ARRAY_TASK_ID;" | cat - case_"$SLURM_ARRAY_TASK_ID".m > temp_"$SLURM_ARRAY_TASK_ID".m && mv temp_"$SLURM_ARRAY_TASK_ID".m case_"$SLURM_ARRAY_TASK_ID".m
echo "clear;" | cat - case_"$SLURM_ARRAY_TASK_ID".m > temp_"$SLURM_ARRAY_TASK_ID".m && mv temp_"$SLURM_ARRAY_TASK_ID".m case_"$SLURM_ARRAY_TASK_ID".m

# Run the sub-job script
matlab -nodisplay < case_"$SLURM_ARRAY_TASK_ID".m

# Delete the sub-job script
rm case_"$SLURM_ARRAY_TASK_ID".m
