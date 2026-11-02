suppressPackageStartupMessages(library(TCGAbiolinks))

options(timeout = 120)

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
task_count <- as.integer(Sys.getenv("GDC_ARRAY_TASK_COUNT", unset = "4"))
if (is.na(task_id) || task_id < 1 || task_id > task_count) {
  stop("SLURM_ARRAY_TASK_ID must be between 1 and ", task_count)
}

analysis_dir <- "/home1/alee4212/490_cluster/analysis_data"
query <- readRDS(file.path(analysis_dir, "brca_methylation_query.rds"))
all_results <- query$results[[1]]
indices <- seq.int(task_id, nrow(all_results), by = task_count)

query_part <- query
query_part$results[[1]] <- all_results[indices, , drop = FALSE]

cat("GDC array task", task_id, "of", task_count,
    "assigned", length(indices), "files\n")
GDCdownload(
  query_part,
  directory = file.path(analysis_dir, "GDCdata"),
  files.per.chunk = 5
)
cat("GDC array task", task_id, "complete\n")
