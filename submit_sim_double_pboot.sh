#!/bin/bash
#SBATCH --job-name=double_pb
#SBATCH --output=job_%A/simulation_%A_%a.out
#SBATCH --error=job_%A/simulation_%A_%a.err
#SBATCH --array=1-50
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=200:00:00
#SBATCH --partition=dscharf-np
#SBATCH --account=dscharf-np

# Load required modules
module load R/4.4

# Create job-specific directory now that we have the job ID
mkdir -p ~/job_${SLURM_ARRAY_JOB_ID}

# Print job info
echo "Starting SLURM job array task ${SLURM_ARRAY_TASK_ID}"
echo "Hostname: $(hostname)"
echo "Date: $(date)"

# Run the R script
Rscript simulation_est_slurm_double_pboot.R ${SLURM_ARRAY_TASK_ID}

echo "Job ${SLURM_ARRAY_TASK_ID} completed at $(date)"
