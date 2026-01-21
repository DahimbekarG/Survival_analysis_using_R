############################################################
# Pan-cancer ZPR1 Kaplan–Meier survival analysis (TCGA)
# Based on LUAD pipeline
# Author: Ganesh Dahimbekar
############################################################

# ----------------------------
# Load libraries
# ----------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(dplyr)
library(tidyr)

# ----------------------------
# TCGA cancer types to analyze
# ----------------------------
tcga_projects <- c(
  "TCGA-LUAD","TCGA-LUSC","TCGA-BRCA","TCGA-COAD",
  "TCGA-READ","TCGA-LIHC","TCGA-STAD","TCGA-PAAD",
  "TCGA-HNSC","TCGA-KIRC"
)

# Create output directory
dir.create("KM_plots", showWarnings = FALSE)

# ----------------------------
# Loop over cancer types
# ----------------------------
for (project in tcga_projects) {
  
  message("\n==============================")
  message("Processing ", project)
  message("==============================")
  
  # ----------------------------
  # 1. Clinical data
  # ----------------------------
  clinical <- tryCatch(
    GDCquery_clinic(project = project, type = "clinical"),
    error = function(e) NULL
  )
  if (is.null(clinical)) {
    message("Clinical data not available for ", project)
    next
  }
  
  clinical <- clinical %>%
    mutate(
      deceased = ifelse(vital_status == "Dead", 1, 0),
      overall_survival = ifelse(
        vital_status == "Dead",
        days_to_death,
        days_to_last_follow_up
      )
    ) %>%
    select(submitter_id, deceased, overall_survival) %>%
    filter(!is.na(overall_survival))
  
  if (nrow(clinical) < 50) {
    message("Too few clinical samples for ", project)
    next
  }
  
  # ----------------------------
  # 2. RNA-seq expression
  # ----------------------------
  query <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = "Primary Tumor"
  )
  
  GDCdownload(query, method = "api")
  se <- GDCprepare(query)
  
  # Use FPKM (as in your LUAD code)
  expr_matrix <- assay(se, "fpkm_unstrand")
  gene_metadata <- as.data.frame(rowData(se))
  
  # ----------------------------
  # 3. Extract ZPR1 expression
  # ----------------------------
  zpr1_id <- gene_metadata %>%
    filter(gene_name == "ZPR1") %>%
    pull(gene_id)
  
  if (length(zpr1_id) == 0) {
    message("ZPR1 not found in ", project)
    next
  }
  
  zpr1_expr <- as.numeric(expr_matrix[zpr1_id, ])
  names(zpr1_expr) <- colnames(expr_matrix)
  
  zpr1_df <- data.frame(
    submitter_id = substr(names(zpr1_expr), 1, 12),
    zpr1_expr = zpr1_expr,
    stringsAsFactors = FALSE
  )
  
  # ----------------------------
  # 4. Merge expression + clinical
  # ----------------------------
  df <- merge(
    zpr1_df,
    clinical,
    by = "submitter_id"
  )
  
  if (nrow(df) < 50) {
    message("Too few merged samples for ", project)
    next
  }
  
  # ----------------------------
  # 5. Median split
  # ----------------------------
  med <- median(df$zpr1_expr, na.rm = TRUE)
  df$group <- ifelse(df$zpr1_expr >= med, "High ZPR1", "Low ZPR1")
  
  # ----------------------------
  # 6. Kaplan–Meier analysis
  # ----------------------------
  fit <- survfit(
    Surv(overall_survival, deceased) ~ group,
    data = df
  )
  
  km_plot <- ggsurvplot(
    fit,
    data = df,
    pval = TRUE,
    risk.table = TRUE,
    conf.int = FALSE,
    xlab = "Time (days)",
    ylab = "Overall Survival Probability",
    title = paste0(project, ": ZPR1 Expression and Overall Survival"),
    legend.title = "ZPR1 Expression",
    legend.labs = c("High ZPR1", "Low ZPR1"),
    risk.table.height = 0.35
  )
  
  # ----------------------------
  # 7. Save plot
  # ----------------------------
  ggsave(
    filename = paste0("KM_plots/", project, "_ZPR1_KM.pdf"),
    plot = km_plot$plot,
    width = 6,
    height = 5
  )
  
  message("KM plot saved for ", project)
}

message("\n===== PAN-CANCER ZPR1 KM ANALYSIS COMPLETE =====")

