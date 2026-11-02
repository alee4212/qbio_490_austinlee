suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(limma)
  library(DESeq2)
  library(data.table)
  library(ggplot2)
})

# Keep individual GDC API attempts short enough to fail fast on a stalled
# connection. GDCdownload is resume-safe and skips files already on disk.
options(timeout = 120)
set.seed(490)

home_dir <- "/home1/alee4212"
analysis_dir <- file.path(home_dir, "490_cluster", "analysis_data")
gdc_dir <- file.path(analysis_dir, "GDCdata")
repo_dir <- file.path(home_dir, "qbio_490_austinlee", "week08_methylation")
output_dir <- file.path(repo_dir, "carc_results")
figure_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(output_dir, "carc_pipeline.log")
zz <- file(log_file, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(zz)
}, add = TRUE)

cat("USC CARC Week 8 methylation pipeline\n")
cat("Started:", format(Sys.time(), tz = "America/Los_Angeles"), "\n")
cat("Host:", Sys.info()[["nodename"]], "\n")
cat("Slurm job:", Sys.getenv("SLURM_JOB_ID", unset = "not set"), "\n")

sample_key <- function(x) substr(x, 1, 15)

cat("\n[1/6] Querying and downloading GDC methylation data\n")
methyl_query_path <- file.path(analysis_dir, "brca_methylation_query.rds")
if (file.exists(methyl_query_path)) {
  methyl_query <- readRDS(methyl_query_path)
} else {
  methyl_query <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "DNA Methylation",
    data.type = "Methylation Beta Value",
    platform = "Illumina Human Methylation 450",
    data.format = "Dataframe"
  )
  saveRDS(methyl_query, methyl_query_path)
}

methyl_se_path <- file.path(analysis_dir, "brca_methylation450_se.rds")
if (file.exists(methyl_se_path)) {
  methylation450 <- readRDS(methyl_se_path)
} else {
  GDCdownload(methyl_query, directory = gdc_dir, method = "client")
  cat("\n[2/6] Preparing the GDC methylation SummarizedExperiment\n")
  methylation450 <- GDCprepare(methyl_query, directory = gdc_dir)
  saveRDS(methylation450, methyl_se_path, compress = FALSE)
}

cat("Methylation object dimensions:",
    paste(dim(methylation450), collapse = " x "), "\n")

cat("\n[3/6] Filtering and modeling methylation\n")
methylation_clinical_all <- as.data.frame(colData(methylation450))
if (!"definition" %in% names(methylation_clinical_all)) {
  stop("The GDC methylation colData does not contain the required definition column.")
}

keep_methyl_samples <- methylation_clinical_all$definition %in%
  c("Primary solid Tumor", "Metastatic")
methylation_clinical <- methylation_clinical_all[keep_methyl_samples, , drop = FALSE]
methylation_clinical$definition <- factor(
  methylation_clinical$definition,
  levels = c("Primary solid Tumor", "Metastatic")
)

betas <- assay(methylation450)[, keep_methyl_samples, drop = FALSE]
stopifnot(identical(colnames(betas), rownames(methylation_clinical)))
cpg_sites_all <- as.data.frame(rowData(methylation450))
if (!"gene" %in% names(cpg_sites_all)) {
  stop("The GDC methylation rowData does not contain the required gene column.")
}

gene_annotation <- as.character(cpg_sites_all$gene)
site_mask <- !is.na(gene_annotation) & gene_annotation != "" &
  !grepl("-", gene_annotation, fixed = TRUE) &
  !grepl(";", gene_annotation, fixed = TRUE) &
  complete.cases(betas)
betas <- betas[site_mask, , drop = FALSE]
cpg_sites <- cpg_sites_all[site_mask, , drop = FALSE]
gene_annotation <- gene_annotation[site_mask]

cat("Methylation sample groups:\n")
print(table(methylation_clinical$definition))
cat("Retained CpG probes:", nrow(betas), "\n")

betas_for_m <- pmin(pmax(betas, 1e-6), 1 - 1e-6)
m_values <- log2(betas_for_m / (1 - betas_for_m))
design <- model.matrix(~ definition, data = methylation_clinical)
methyl_fit <- eBayes(lmFit(m_values, design))
coefficient_name <- "definitionMetastatic"

methylation_results <- data.table(
  probe_id = rownames(betas),
  gene_name = gene_annotation,
  m_value_difference = methyl_fit$coefficients[, coefficient_name],
  p_value = methyl_fit$p.value[, coefficient_name]
)
methylation_results[, adj_p_value := p.adjust(p_value, method = "BY")]
methylation_results[, neg_log10_adj_p :=
  -log10(pmax(adj_p_value, .Machine$double.xmin))]
methylation_results[, significant :=
  adj_p_value < 0.05 & abs(m_value_difference) >= 1]
setorder(methylation_results, adj_p_value)
fwrite(methylation_results,
       file.path(output_dir, "methylation_results_carc.csv.gz"),
       compress = "gzip")
fwrite(head(methylation_results, 100),
       file.path(output_dir, "methylation_top100_carc.csv"))

volcano <- ggplot(
  methylation_results,
  aes(m_value_difference, neg_log10_adj_p, color = significant)
) +
  geom_point(alpha = 0.35, size = 0.55) +
  geom_vline(xintercept = c(-1, 1), color = "#8B0000", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), color = "#8B0000", linetype = "dashed") +
  scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "#2C7FB8")) +
  labs(
    title = "TCGA-BRCA Differential Methylation (USC CARC)",
    subtitle = "Metastatic versus primary solid tumor",
    x = "M-value difference (metastatic - primary)",
    y = "-log10(BY-adjusted p-value)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
ggsave(file.path(figure_dir, "01_methylation_volcano_carc.png"),
       volcano, width = 8, height = 6, dpi = 300)

rm(m_values, betas_for_m, methyl_fit, methylation450)
gc()

cat("\n[4/6] Preparing the existing CARC GDC RNA download\n")
rna_query_path <- file.path(analysis_dir, "brca_rna_query.rds")
if (file.exists(rna_query_path)) {
  rna_query <- readRDS(rna_query_path)
} else {
  rna_query <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  saveRDS(rna_query, rna_query_path)
}

rna_se_path <- file.path(analysis_dir, "brca_rna_star_counts_se.rds")
if (file.exists(rna_se_path)) {
  rna_se <- readRDS(rna_se_path)
} else {
  rna_se <- GDCprepare(rna_query, directory = gdc_dir)
  saveRDS(rna_se, rna_se_path, compress = FALSE)
}

assay_to_use <- if ("unstranded" %in% assayNames(rna_se))
  "unstranded" else assayNames(rna_se)[[1]]
rna_counts_all <- assay(rna_se, assay_to_use)
rna_gene_data <- as.data.frame(rowData(rna_se))
rna_gene_names <- if ("gene_name" %in% names(rna_gene_data))
  as.character(rna_gene_data$gene_name) else rownames(rna_gene_data)

# Keep one aliquot per biological sample and exactly match the 15-character
# TCGA participant + sample-type key to the methylation clinical data.
rna_keys <- sample_key(colnames(rna_counts_all))
unique_rna <- !duplicated(rna_keys)
rna_counts_all <- rna_counts_all[, unique_rna, drop = FALSE]
rna_keys <- rna_keys[unique_rna]
methyl_keys <- sample_key(rownames(methylation_clinical))
matched_methyl <- match(rna_keys, methyl_keys)
matched <- !is.na(matched_methyl)
rna_counts <- rna_counts_all[, matched, drop = FALSE]
matched_methyl <- matched_methyl[matched]

rna_clinical <- data.frame(
  definition = factor(
    methylation_clinical$definition[matched_methyl],
    levels = c("Primary solid Tumor", "Metastatic")
  ),
  row.names = colnames(rna_counts)
)
stopifnot(identical(colnames(rna_counts), rownames(rna_clinical)))

keep_rna_genes <- rowSums(rna_counts) >= 20
rna_counts <- rna_counts[keep_rna_genes, , drop = FALSE]
rna_gene_names <- rna_gene_names[keep_rna_genes]
storage.mode(rna_counts) <- "integer"

cat("RNA sample groups:\n")
print(table(rna_clinical$definition))
cat("Genes after total-count filter:", nrow(rna_counts), "\n")

cat("\n[5/6] Running DESeq2 and intersecting candidate genes\n")
dds <- DESeqDataSetFromMatrix(
  countData = rna_counts,
  colData = rna_clinical,
  design = ~ definition
)
dds <- DESeq(dds, quiet = TRUE)
de_results <- as.data.table(as.data.frame(results(
  dds,
  contrast = c("definition", "Metastatic", "Primary solid Tumor")
)), keep.rownames = "gene_id")
de_results[, gene_name := rna_gene_names[match(gene_id, rownames(rna_counts))]]
setcolorder(de_results,
            c("gene_id", "gene_name",
              setdiff(names(de_results), c("gene_id", "gene_name"))))
de_results <- de_results[order(is.na(padj), padj)]
fwrite(de_results, file.path(output_dir, "deseq2_results_carc.csv.gz"),
       compress = "gzip")
fwrite(head(de_results[!is.na(padj)], 100),
       file.path(output_dir, "deseq2_top100_carc.csv"))

downregulated <- unique(de_results[
  !is.na(padj) & padj < 0.05 & log2FoldChange < -1 & !is.na(gene_name),
  gene_name
])
hypermethylated <- unique(methylation_results[
  m_value_difference > 1.5 & !is.na(gene_name), gene_name
])
interest_genes <- intersect(downregulated, hypermethylated)
if (length(interest_genes) == 0) {
  stop("No genes passed both the RNA and methylation thresholds in the CARC data.")
}

candidate_rna <- de_results[gene_name %in% interest_genes]
candidate_rna <- candidate_rna[order(padj, -abs(log2FoldChange))]
selected_gene <- candidate_rna$gene_name[[1]]
selected_gene_id <- candidate_rna[gene_name == selected_gene, gene_id][[1]]
cat("Interest genes:", paste(interest_genes, collapse = ", "), "\n")
cat("Selected gene:", selected_gene, "\n")

normalized_counts <- counts(dds, normalized = TRUE)[selected_gene_id, ]
expression_plot_data <- data.table(
  sample = names(normalized_counts),
  definition = rna_clinical[names(normalized_counts), "definition"],
  normalized_count = as.numeric(normalized_counts)
)
fwrite(expression_plot_data,
       file.path(output_dir, "selected_gene_expression_plot_data_carc.csv"))

expression_plot <- ggplot(
  expression_plot_data,
  aes(definition, log10(normalized_count + 1), fill = definition)
) +
  geom_boxplot(width = 0.6, outlier.alpha = 0.3) +
  geom_jitter(width = 0.12, alpha = 0.3, size = 0.8) +
  scale_fill_manual(values = c(
    "Primary solid Tumor" = "#377EB8", "Metastatic" = "#E41A1C"
  )) +
  labs(
    title = paste0(selected_gene, " RNA Expression (USC CARC)"),
    subtitle = "TCGA-BRCA metastatic versus primary solid tumor",
    x = NULL,
    y = "log10(DESeq2 normalized count + 1)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
ggsave(file.path(figure_dir, "02_selected_gene_rna_expression_carc.png"),
       expression_plot, width = 7, height = 5, dpi = 300)

selected_probe_rows <- which(gene_annotation == selected_gene)
beta_means <- data.table(
  probe_id = rownames(betas)[selected_probe_rows],
  primary_mean_beta = rowMeans(
    betas[selected_probe_rows,
          methylation_clinical$definition == "Primary solid Tumor", drop = FALSE]
  ),
  metastatic_mean_beta = rowMeans(
    betas[selected_probe_rows,
          methylation_clinical$definition == "Metastatic", drop = FALSE]
  )
)
beta_means[, probe_order := seq_len(.N)]
fwrite(beta_means,
       file.path(output_dir, "selected_gene_methylation_plot_data_carc.csv"))

beta_long <- melt(
  beta_means,
  id.vars = c("probe_id", "probe_order"),
  measure.vars = c("primary_mean_beta", "metastatic_mean_beta"),
  variable.name = "definition",
  value.name = "mean_beta"
)
beta_long[, definition := factor(
  definition,
  levels = c("primary_mean_beta", "metastatic_mean_beta"),
  labels = c("Primary solid Tumor", "Metastatic")
)]

methylation_gene_plot <- ggplot(
  beta_long,
  aes(probe_order, mean_beta, fill = definition)
) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  scale_fill_manual(values = c(
    "Primary solid Tumor" = "#377EB8", "Metastatic" = "#E41A1C"
  )) +
  labs(
    title = paste0(selected_gene, " CpG Methylation (USC CARC)"),
    subtitle = "Mean beta value at each retained CpG probe",
    x = "CpG probe (ordered)",
    y = "Mean beta value",
    fill = "Sample type"
  ) +
  theme_bw(base_size = 12)
ggsave(file.path(figure_dir, "03_selected_gene_cpg_methylation_carc.png"),
       methylation_gene_plot, width = 9, height = 5, dpi = 300)

cat("\n[6/6] Writing CARC summary and provenance\n")
summary_table <- data.table(
  metric = c(
    "Execution platform", "CARC host", "Slurm job ID",
    "Methylation primary samples", "Methylation metastatic samples",
    "Retained CpG probes", "Significant CpG probes",
    "RNA primary samples", "RNA metastatic samples",
    "Genes after count filter", "Significant downregulated genes",
    "Effect-size hypermethylated genes", "Intersecting genes", "Selected gene"
  ),
  value = c(
    "USC CARC Discovery", Sys.info()[["nodename"]],
    Sys.getenv("SLURM_JOB_ID", unset = "not set"),
    sum(methylation_clinical$definition == "Primary solid Tumor"),
    sum(methylation_clinical$definition == "Metastatic"),
    nrow(betas), sum(methylation_results$significant),
    sum(rna_clinical$definition == "Primary solid Tumor"),
    sum(rna_clinical$definition == "Metastatic"),
    nrow(rna_counts), length(downregulated), length(hypermethylated),
    length(interest_genes), selected_gene
  )
)
fwrite(summary_table, file.path(output_dir, "analysis_summary_carc.csv"))
fwrite(data.table(gene_name = interest_genes),
       file.path(output_dir, "interest_genes_carc.csv"))
writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "sessionInfo_carc.txt"))

cat("Completed:", format(Sys.time(), tz = "America/Los_Angeles"), "\n")
cat("CARC analysis complete.\n")
