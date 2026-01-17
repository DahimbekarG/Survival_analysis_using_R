#!/usr/bin/env Rscript
# ============================================================
# Pan-cancer ZPR1 survival analysis using TCGA RNA-seq data
# Author: Ganesh Dahimbekar
# Method: DESeq2 VST normalization + KM + Cox models
# Robustized: clinical fallback, deduplication, non-count assay fallback
# ============================================================

# ---- Helpers to install packages non-interactively ----
install_if_missing_cran <- function(pkgs){
  to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if(length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org", dependencies = TRUE)
}
install_if_missing_bioc <- function(pkgs){
  if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if(length(to_install)) BiocManager::install(to_install, ask = FALSE, update = FALSE)
}

cran_pkgs <- c("dplyr","tidyr","ggplot2","survival","survminer")
bioc_pkgs <- c("TCGAbiolinks","SummarizedExperiment","DESeq2")
install_if_missing_cran(cran_pkgs)
install_if_missing_bioc(bioc_pkgs)

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(survival)
  library(survminer)
  library(ggplot2)
})

message("===== TCGA PAN-CANCER ZPR1 SURVIVAL ANALYSIS =====")

# Projects to analyze (edit as desired)
tcga_projects <- c(
  "TCGA-LUAD", "TCGA-LUSC", "TCGA-BRCA",
  "TCGA-COAD", "TCGA-READ", "TCGA-LIHC",
  "TCGA-STAD", "TCGA-PAAD", "TCGA-HNSC",
  "TCGA-KIRC"
)

# Core function per project
run_zpr1_survival <- function(project, min_samples = 30){
  message(sprintf("\n--- Processing %s ---", project))

  # 1) Try clinic via API (non-fatal)
  clin_raw <- tryCatch(GDCquery_clinic(project = project, type = "clinical"), error = function(e) { message("GDCquery_clinic failed: ", e$message); NULL })

  # 2) Query expression
  query <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = "Primary Tumor"
  )

  GDCdownload(query, method = "api", files.per.chunk = 20)
  se <- GDCprepare(query)

  # 3) Build clinical table: prefer clinic_raw, else colData fallback
  clinical <- NULL
  if (!is.null(clin_raw)) {
    clinical <- as.data.frame(clin_raw, stringsAsFactors = FALSE)
    # try to detect id column and survival cols
    id_candidates <- c("submitter_id","bcr_patient_barcode","patient_id","case_id")
    id_col <- intersect(id_candidates, colnames(clinical))
    if (length(id_col) > 0) id_col <- id_col[1] else stop("No patient id column in clinical data")
    clinical <- clinical %>%
      mutate(deceased = ifelse(tolower(vital_status) == "dead", 1L, 0L),
             overall_survival = as.numeric(ifelse(deceased == 1L, days_to_death, days_to_last_follow_up))) %>%
      dplyr::rename(submitter_id = !!rlang::sym(id_col)) %>%
      dplyr::select(submitter_id, deceased, overall_survival) %>%
      filter(!is.na(overall_survival) & overall_survival > 0)
  } else {
    message("Using colData(se) to construct clinical table")
    cd <- as.data.frame(colData(se), stringsAsFactors = FALSE)
    # detect ids and survival columns
    id_candidates_cd <- c("submitter_id","bcr_patient_barcode","patient","patient_id","case_id","barcode","sample_submitter_id")
    id_col_cd <- intersect(id_candidates_cd, colnames(cd))
    if (length(id_col_cd) == 0) stop("No identifier column found in colData(se)")
    id_col_cd <- id_col_cd[1]

    vs_col <- intersect(c("vital_status","vital.status"), colnames(cd))
    dtd_col <- intersect(c("days_to_death"), colnames(cd))
    dtlfu_col <- intersect(c("days_to_last_follow_up"), colnames(cd))
    if (length(vs_col) == 0) stop("No vital_status in colData(se)")
    vs_col <- vs_col[1]
    dtd_col <- if(length(dtd_col)>0) dtd_col[1] else NA
    dtlfu_col <- if(length(dtlfu_col)>0) dtlfu_col[1] else NA

    deceased_vec <- ifelse(tolower(as.character(cd[[vs_col]])) == "dead", 1L, 0L)
    overall_surv_vec <- rep(NA_real_, nrow(cd))
    if (!is.na(dtd_col)) overall_surv_vec[!is.na(cd[[dtd_col]])] <- as.numeric(cd[[dtd_col]][!is.na(cd[[dtd_col]])])
    if (!is.na(dtlfu_col)) {
      sel <- is.na(overall_surv_vec) & !is.na(cd[[dtlfu_col]])
      overall_surv_vec[sel] <- as.numeric(cd[[dtlfu_col]][sel])
    }

    clinical <- data.frame(submitter_id = as.character(cd[[id_col_cd]]),
                           deceased = deceased_vec, overall_survival = overall_surv_vec, stringsAsFactors = FALSE)
    clinical <- clinical[!is.na(clinical$overall_survival) & clinical$overall_survival>0, , drop=FALSE]
  }

  if (nrow(clinical) < min_samples) { message("Not enough clinical samples: ", nrow(clinical)); return(NULL) }

  # 4) Expression assay selection (prefer counts-like)
  assays_available <- assayNames(se)
  pref <- intersect(c("unstranded","raw_count","HTSeq - Counts","counts","STAR - Counts"), assays_available)
  assay_name <- if (length(pref)>0) pref[1] else assays_available[1]
  message("Using assay: ", assay_name)
  expr_mat <- assay(se, assay_name)

  gene_md <- as.data.frame(rowData(se))
  # standardize gene name column
  if (!"gene_name" %in% colnames(gene_md)) {
    alt <- intersect(c("external_gene_name","gene_symbol","symbol"), colnames(gene_md))
    if (length(alt)>0) colnames(gene_md)[which(colnames(gene_md)==alt[1])] <- "gene_name"
  }
  if (!"gene_name" %in% colnames(gene_md)) stop("No gene name column in rowData")

  zpr1_rows <- which(toupper(gene_md$gene_name) == "ZPR1")
  if (length(zpr1_rows) == 0) { message("ZPR1 not found in project: ", project); return(NULL) }

  # If multiple aliquots per patient, aggregate to patient-level first
  col_names <- colnames(expr_mat)
  patient_id <- substr(col_names, 1, 12)
  if (any(duplicated(patient_id))) {
    message("Multiple aliquots per patient detected; aggregating expression by patient (mean) to avoid many-to-many joins")
    # Aggregate all genes to patient-level (mean) to preserve one sample per patient
    agg_mat <- t(sapply(split(seq_len(ncol(expr_mat)), patient_id), function(cols) rowMeans(expr_mat[, cols, drop=FALSE], na.rm = TRUE)))
    # agg_mat rows are patients, columns are genes -> transpose back
    expr_mat_pat <- t(agg_mat)
    colnames(expr_mat_pat) <- names(split(seq_len(ncol(expr_mat)), patient_id))
  } else {
    expr_mat_pat <- expr_mat
  }

  # locate zpr1 row in patient-level matrix
  if (length(zpr1_rows) > 1) {
    # sum counts if counts-like; else mean
    if (all(expr_mat == floor(expr_mat), na.rm = TRUE)) {
      zpr1_vals <- colSums(expr_mat[zpr1_rows, , drop=FALSE], na.rm = TRUE)
    } else {
      zpr1_vals <- colMeans(expr_mat[zpr1_rows, , drop=FALSE], na.rm = TRUE)
    }
    # if we aggregated to patients, re-aggregate zpr1_vals accordingly
    zpr1_patient <- tapply(zpr1_vals, patient_id, mean)
  } else {
    zpr1_vals <- expr_mat[zpr1_rows, ]
    zpr1_patient <- tapply(as.numeric(zpr1_vals), patient_id, mean)
  }

  # If counts-like, run DESeq2 VST on patient-level counts; else use log2(TPM/FPKM + 1)
  is_counts_like <- all(expr_mat == floor(expr_mat), na.rm = TRUE)
  if (is_counts_like) {
    message("Counts-like assay detected. Running DESeq2 VST on patient-level matrix (this may be slow).")
    # Build patient-level counts matrix by summing per patient if duplicates existed originally
    if (any(duplicated(patient_id))) {
      # sum counts per gene for each patient
      counts_pat <- t(sapply(split(seq_len(ncol(expr_mat)), patient_id), function(cols) rowSums(expr_mat[, cols, drop=FALSE], na.rm = TRUE)))
      counts_pat <- t(counts_pat)
      colnames(counts_pat) <- names(split(seq_len(ncol(expr_mat)), patient_id))
    } else {
      counts_pat <- expr_mat
    }
    dds <- DESeqDataSetFromMatrix(countData = counts_pat, colData = DataFrame(row.names = colnames(counts_pat)), design = ~1)
    dds <- estimateSizeFactors(dds)
    vst_mat <- vst(dds, blind = TRUE)
    # normalized zpr1 expression (patient-level)
    zpr1_norm <- as.numeric(assay(vst_mat)[which(rownames(assay(vst_mat)) %in% rownames(expr_mat)[zpr1_rows])[1], ])
    names(zpr1_norm) <- colnames(assay(vst_mat))
  } else {
    message("Non-count assay detected; using log2(x+1) transform for ZPR1")
    zpr1_norm <- log2(zpr1_patient + 1)
    # names are patient ids
  }

  # Build expression data.frame and merge with clinical
  expr_df <- data.frame(submitter_id = names(zpr1_norm), zpr1_expr = as.numeric(zpr1_norm), stringsAsFactors = FALSE)
  df <- left_join(expr_df, clinical, by = "submitter_id") %>% filter(!is.na(overall_survival))

  if (nrow(df) < min_samples) { message("Not enough merged samples: ", nrow(df)); return(NULL) }

  # Median split
  med <- median(df$zpr1_expr, na.rm = TRUE)
  df$group <- factor(ifelse(df$zpr1_expr >= med, "High ZPR1", "Low ZPR1"), levels = c("High ZPR1","Low ZPR1"))

  # Cox models
  cox_cont <- coxph(Surv(overall_survival, deceased) ~ zpr1_expr, data = df)
  cox_bin  <- coxph(Surv(overall_survival, deceased) ~ group, data = df)

  res <- data.frame(
    project = project,
    n = nrow(df),
    HR_cont = exp(coef(cox_cont)),
    p_cont = summary(cox_cont)$coefficients[,
                                            "Pr(>|z|)"][[1]],
    HR_bin = exp(coef(cox_bin)),
    p_bin = summary(cox_bin)$coefficients[,
                                         "Pr(>|z|)"][[1]],
    stringsAsFactors = FALSE
  )
  return(res)
}

# Run across projects
results_list <- lapply(tcga_projects, function(p) tryCatch(run_zpr1_survival(p), error = function(e) { message("Error in project ", p, ": ", e$message); NULL }))
pan_results <- bind_rows(results_list)
write.csv(pan_results, "ZPR1_pan_cancer_survival_results.csv", row.names = FALSE)
print(pan_results)

# Forest plot if any results
if (nrow(pan_results) > 0) {
  forest <- ggplot(pan_results, aes(x = reorder(project, HR_cont), y = HR_cont)) +
    geom_point(size = 3) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    coord_flip() +
    ylab("Hazard Ratio (ZPR1 expression)") + xlab("Cancer type") + theme_bw() + ggtitle("Pan-cancer ZPR1 survival association (TCGA)")
  ggsave("ZPR1_pan_cancer_forest_plot.pdf", forest, width = 7, height = 5)
}

message("===== ANALYSIS COMPLETE =====")
