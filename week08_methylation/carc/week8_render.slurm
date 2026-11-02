#!/bin/bash
#SBATCH --job-name=week8_render
#SBATCH --account=rohs_1070
#SBATCH --partition=main
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/home1/alee4212/qbio_490_austinlee/week08_methylation/carc/week8_render_%j.out
#SBATCH --error=/home1/alee4212/qbio_490_austinlee/week08_methylation/carc/week8_render_%j.err

set -euo pipefail

module purge
module load ver/2506 gcc/14.3.0 r/4.5.0 pandoc/3.9.0.2

export R_LIBS_USER=/home1/alee4212/R/x86_64-pc-linux-gnu-library/4.5

Rscript -e 'if (!requireNamespace("rmarkdown", quietly = TRUE)) install.packages("rmarkdown", repos = "https://cloud.r-project.org")'

cd /home1/alee4212/qbio_490_austinlee/week08_methylation
Rscript -e 'rmarkdown::render("Methylation_Homework.Rmd", output_file = "Methylation_Homework.nb.html", quiet = FALSE)'
