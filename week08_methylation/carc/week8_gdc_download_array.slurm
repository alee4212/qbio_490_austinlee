#!/bin/bash
#SBATCH --job-name=week8_gdc
#SBATCH --account=rohs_1070
#SBATCH --partition=main
#SBATCH --array=1-4
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/home1/alee4212/qbio_490_austinlee/week08_methylation/carc/week8_gdc_%A_%a.out
#SBATCH --error=/home1/alee4212/qbio_490_austinlee/week08_methylation/carc/week8_gdc_%A_%a.err

set -euo pipefail

module purge
module load ver/2506 gcc/14.3.0 r/4.5.0

export R_LIBS_USER=/home1/alee4212/R/x86_64-pc-linux-gnu-library/4.5
export GDC_ARRAY_TASK_COUNT=4

task_tmp=/home1/alee4212/490_cluster/gdc_parts/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$task_tmp"
cd "$task_tmp"

Rscript /home1/alee4212/qbio_490_austinlee/week08_methylation/carc/week8_gdc_download_part.R
