############################################################
# Pan-cancer ZPR1 Kaplan–Meier survival analysis (TCGA)
# Based on LUAD pipeline
# Author: Ganesh Dahimbekar
############################################################

# 1. Load Essential Libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(dplyr)
library(httr)

# 2. Setup Working Environment
# Critical for Windows: Use a very short path to avoid GDCprepare errors
work_dir <- "C:/TCGA_ZPR1"
if(!dir.exists(work_dir)) dir.create(work_dir)
setwd(work_dir)
dir.create("KM_plots", showWarnings = FALSE)

# Increase timeout for large 2026 GDC files
options(timeout = 600) 

tcga_projects <- c("TCGA-LUAD", "TCGA-LUSC", "TCGA-BRCA", "TCGA-COAD", 
                   "TCGA-READ", "TCGA-LIHC", "TCGA-STAD", "TCGA-PAAD", 
                   "TCGA-HNSC", "TCGA-KIRC")

# 3. Processing Loop
for (project in tcga_projects) {
  message("\nProcessing: ", project)
  
  # --- Step A: Clinical Data ---
  clinical <- tryCatch({ GDCquery_clinic(project = project, type = "clinical") }, 
                       error = function(e) return(NULL))
  
  if (is.null(clinical) || nrow(clinical) == 0) next
  
  clinical_clean <- clinical %>%
    mutate(deceased = ifelse(vital_status == "Dead", 1, 0),
           overall_survival = ifelse(vital_status == "Dead", days_to_death, days_to_last_follow_up)) %>%
    filter(!is.na(overall_survival) & overall_survival > 0) %>%
    distinct(submitter_id, .keep_all = TRUE) # Prevents duplication errors
  
  # --- Step B: Gene Expression Query ---
  query <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = "Primary Tumor"
  )
  
  # --- Step C: Robust Download (Fixes "Truncated tar" error) ---
  # Retries with small chunks or the official GDC-Client if API fails
  tryCatch({
    GDCdownload(query, method = "api", files.per.chunk = 15)
  }, error = function(e) {
    message("API failed for ", project, ". Attempting 'client' method...")
    GDCdownload(query, method = "client")
  })
  
  # --- Step D: Prepare & Extract ZPR1 ---
  # Explicitly naming the directory avoids the GDCprepare "not found" error
  se <- GDCprepare(query, directory = "GDCdata")
  
  # Dynamic assay detection (FPKM name may vary slightly in 2026)
  assay_target <- grep("fpkm_unstrand", assayNames(se), value = TRUE)[1]
  if(is.na(assay_target)) assay_target <- "fpkm_unstranded"
  
  expr_matrix <- assay(se, assay_target)
  zpr1_id <- rownames(rowData(se)[rowData(se)$gene_name == "ZPR1", , drop=FALSE])
  
  if (length(zpr1_id) == 0) { message("ZPR1 missing in ", project); next }
  
  zpr1_df <- data.frame(
    submitter_id = substr(colnames(expr_matrix), 1, 12),
    zpr1_expr = as.numeric(expr_matrix[zpr1_id, ]),
    stringsAsFactors = FALSE
  ) %>% group_by(submitter_id) %>% summarize(zpr1_expr = mean(zpr1_expr))
  
  # --- Step E: Merge & Survival Analysis ---
  df_final <- inner_join(zpr1_df, clinical_clean, by = "submitter_id")
  if(nrow(df_final) < 20) next
  
  df_final$group <- ifelse(df_final$zpr1_expr >= median(df_final$zpr1_expr), "High", "Low")
  fit <- survfit(Surv(overall_survival, deceased) ~ group, data = df_final)
  
  # --- Step F: Save Results ---
  p <- ggsurvplot(fit, data = df_final, pval = TRUE, risk.table = TRUE,
                  title = paste(project, "ZPR1 Survival"))
  
  pdf(file.path("KM_plots", paste0(project, "_ZPR1.pdf")), onefile = FALSE)
  print(p)
  dev.off()
  
  # --- Step G: Memory Cleanup (Prevents R session crashes) ---
  rm(se, expr_matrix, query, df_final, clinical, zpr1_df)
  gc()
}

message("\nAnalysis Complete. Plots saved to: ", getwd(), "/KM_plots")

message("\n===== PAN-CANCER ZPR1 KM ANALYSIS COMPLETE =====")

