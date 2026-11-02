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
module load ver/2506 gcc/14.3.0 r/4.5.0 gdc-client/2.3

export R_LIBS_USER=/home1/alee4212/R/x86_64-pc-linux-gnu-library/4.5

Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org"); BiocManager::install(c("limma", "DESeq2"), ask=FALSE, update=FALSE)'

pipeline=/home1/alee4212/qbio_490_austinlee/week08_methylation/carc/week8_carc_pipeline.R
pipeline_log=/home1/alee4212/qbio_490_austinlee/week08_methylation/carc_results/carc_pipeline.log

for attempt in 1 2 3 4 5 6 7 8; do
  echo "CARC pipeline attempt ${attempt}/8"
  if Rscript "$pipeline"; then
    exit 0
  fi

  if grep -q "GDC server down\|api.gdc.cancer.gov/status.*Timeout" "$pipeline_log"; then
    echo "Temporary GDC outage detected; retrying in 60 seconds."
    sleep 60
  else
    echo "Non-network pipeline failure; not retrying."
    exit 1
  fi
done

echo "GDC remained unavailable after eight attempts."
exit 1
