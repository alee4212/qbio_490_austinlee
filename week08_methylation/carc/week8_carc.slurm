#!/bin/bash
#SBATCH --job-name=week8_methylation
#SBATCH --account=rohs_1070
#SBATCH --partition=main
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --output=/home1/alee4212/qbio_490_austinlee/week08_methylation/carc/week8_carc_%j.out
#SBATCH --error=/home1/alee4212/qbio_490_austinlee/week08_methylation/carc/week8_carc_%j.err

set -euo pipefail

module purge
module load ver/2506 gcc/14.3.0 r/4.5.0

export R_LIBS_USER=/home1/alee4212/R/x86_64-pc-linux-gnu-library/4.5

Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org"); BiocManager::install(c("limma", "DESeq2"), ask=FALSE, update=FALSE)'

Rscript /home1/alee4212/qbio_490_austinlee/week08_methylation/carc/week8_carc_pipeline.R
