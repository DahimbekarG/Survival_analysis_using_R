#!/usr/bin/env Rscript
# ============================================================
# Pan-cancer ZPR1 survival analysis using TCGA RNA-seq data
# Author: Ganesh Dahimbekar
# Method: DESeq2 VST normalization + KM + Cox models
# ============================================================

# ---- Helpers to install packages ----
install_if_missing_cran <- function(pkgs){
  to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if(length(to_install))
    install.packages(to_install, repos = "https://cloud.r-project.org", dependencies = TRUE)
}
install_if_missing_bioc <- function(pkgs){
  if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if(length(to_install))
    BiocManager::install(to_install, ask = FALSE, update = FALSE)
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

dir.create("KM_plots", showWarnings = FALSE)

message("===== TCGA PAN-CANCER ZPR1 SURVIVAL ANALYSIS =====")

# Cancer projects
tcga_projects <- c(
  "TCGA-LUAD","TCGA-LUSC","TCGA-BRCA","TCGA-COAD",
  "TCGA-READ","TCGA-LIHC","TCGA-STAD","TCGA-PAAD",
  "TCGA-HNSC","TCGA-KIRC"
)

# ------------------------------------------------------------
# Core function
# ------------------------------------------------------------
run_zpr1_survival <- function(project, min_samples = 30){

  message("\n--- Processing ", project, " ---")

  # ---- Clinical ----
  clinical <- tryCatch(
    GDCquery_clinic(project = project, type = "clinical"),
    error = function(e) NULL
  )
  if (is.null(clinical)) return(NULL)

  clinical <- as.data.frame(clinical) %>%
    mutate(
      deceased = ifelse(tolower(vital_status) == "dead", 1L, 0L),
      overall_survival = ifelse(
        deceased == 1L,
        days_to_death,
        days_to_last_follow_up
      )
    ) %>%
    select(submitter_id, deceased, overall_survival) %>%
    filter(!is.na(overall_survival) & overall_survival > 0)

  if (nrow(clinical) < min_samples) return(NULL)

  # ---- Expression ----
  query <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = "Primary Tumor"
  )
  GDCdownload(query, method = "api", files.per.chunk = 20)
  se <- GDCprepare(query)

  counts <- assay(se, "unstranded")
  gene_md <- as.data.frame(rowData(se))
  if (!"gene_name" %in% colnames(gene_md)) return(NULL)

  zpr1_rows <- which(gene_md$gene_name == "ZPR1")
  if (length(zpr1_rows) == 0) return(NULL)

  # ---- Aggregate to patient level ----
  patient_id <- substr(colnames(counts), 1, 12)
  counts_pat <- t(sapply(split(seq_len(ncol(counts)), patient_id),
                         function(i) rowSums(counts[, i, drop=FALSE])))
  counts_pat <- t(counts_pat)

  # ---- DESeq2 VST ----
  dds <- DESeqDataSetFromMatrix(
    countData = counts_pat,
    colData = DataFrame(row.names = colnames(counts_pat)),
    design = ~1
  )
  dds <- estimateSizeFactors(dds)
  vst_mat <- vst(dds, blind = TRUE)

  zpr1_expr <- assay(vst_mat)[zpr1_rows[1], ]
  expr_df <- data.frame(
    submitter_id = colnames(vst_mat),
    zpr1_expr = as.numeric(zpr1_expr)
  )

  # ---- Merge ----
  df <- left_join(expr_df, clinical, by = "submitter_id")
  if (nrow(df) < min_samples) return(NULL)

  # ---- Median split ----
  med <- median(df$zpr1_expr, na.rm = TRUE)
  df$group <- factor(
    ifelse(df$zpr1_expr >= med, "High ZPR1", "Low ZPR1"),
    levels = c("High ZPR1","Low ZPR1")
  )

  # ======================
  # Kaplan–Meier plot
  # ======================
  fit <- survfit(Surv(overall_survival, deceased) ~ group, data = df)

  km_plot <- ggsurvplot(
    fit,
    data = df,
    pval = TRUE,
    risk.table = TRUE,
    conf.int = FALSE,
    xlab = "Time (days)",
    ylab = "Overall survival probability",
    title = paste0("ZPR1 and Overall Survival in ", project),
    legend.title = "Expression",
    legend.labs = c("High ZPR1","Low ZPR1"),
    risk.table.height = 0.25
  )

  ggsave(
    filename = file.path("KM_plots", paste0("KM_ZPR1_", project, ".pdf")),
    plot = km_plot$plot,
    width = 7, height = 6
  )

  # ---- Cox ----
  cox_cont <- coxph(Surv(overall_survival, deceased) ~ zpr1_expr, data = df)
  cox_bin  <- coxph(Surv(overall_survival, deceased) ~ group, data = df)

  data.frame(
    project = project,
    n = nrow(df),
    HR_cont = exp(coef(cox_cont)),
    p_cont = summary(cox_cont)$coefficients[,"Pr(>|z|)"],
    HR_bin = exp(coef(cox_bin)),
    p_bin = summary(cox_bin)$coefficients[,"Pr(>|z|)"]
  )
}

# ------------------------------------------------------------
# Run pan-cancer
# ------------------------------------------------------------
results <- lapply(tcga_projects, function(p)
  tryCatch(run_zpr1_survival(p), error = function(e) NULL))

pan_results <- bind_rows(results)
write.csv(pan_results, "ZPR1_pan_cancer_survival_results.csv", row.names = FALSE)
print(pan_results)

# ------------------------------------------------------------
# Forest plot
# ------------------------------------------------------------
if (nrow(pan_results) > 0) {
  forest <- ggplot(pan_results, aes(x = reorder(project, HR_cont), y = HR_cont)) +
    geom_point(size = 3) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    coord_flip() +
    ylab("Hazard Ratio (ZPR1 expression)") +
    xlab("Cancer type") +
    theme_bw() +
    ggtitle("Pan-cancer ZPR1 survival association (TCGA)")

  ggsave("ZPR1_pan_cancer_forest_plot.pdf", forest, width = 7, height = 5)
}

message("===== ANALYSIS COMPLETE =====")
