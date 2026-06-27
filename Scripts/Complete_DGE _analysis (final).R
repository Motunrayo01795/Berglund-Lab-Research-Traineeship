# ==============================================================================
# SCRIPT: Myotonic Dystrophy (DM1) DGE, Rescue, and Off-Target Toxicity Pipeline
# DESCRIPTION: Comprehensive DESeq2 analysis evaluating the therapeutic efficacy 
#              and off-target transcriptomic effects of ASO, SAHA, and Combo 
#              treatments in DM1 patient-derived myoblasts.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. SETUP & LIBRARIES
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  # Core Data Manipulation
  library(tidyverse)
  
  # Differential Expression
  library(DESeq2)
  
  # Annotation & Pathway Analysis
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(clusterProfiler)
  
  # Visualization
  library(ggplot2)
  library(ggvenn)
  library(UpSetR)
  library(grid)
})

# ------------------------------------------------------------------------------
# 2. CONFIGURATION & DIRECTORIES
# ------------------------------------------------------------------------------
# Set input and output directories
FASTQ_DIR <- "Z:\\common\\Sequencing_Data\\projects\\Hoque_Myoblast_01212025\\For_NAS\\FASTQ\\Trimmed_reads"
OUT_DIR   <- "C:\\Users\\USER\\Documents\\Berglund\\Phase_two\\Phase 2\\Outputs"

# Ensure output directory exists
if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# Define Sample Groups
samples_ctrl  <- c("CTRL2_MB_DMSO_A_S1", "CTRL2_MB_DMSO_B_S2", "CTRL2_MB_DMSO_C_S3")
samples_dm1   <- c("DM05_MB_DMSO_A_S4", "DM05_MB_DMSO_B_S5", "DM05_MB_DMSO_C_S6")
samples_saha  <- c("DM05_MB_SAHA_A_S7", "DM05_MB_SAHA_B_S8", "DM05_MB_SAHA_C_S9")
samples_aso   <- c("DM05_MB_ASO_A_S10", "DM05_MB_ASO_B_S11", "DM05_MB_ASO_C_S12")
samples_combo <- c("DM05_MB_ASO_SAHA_A_S13", "DM05_MB_ASO_SAHA_B_S14", "DM05_MB_ASO_SAHA_C_S15")

all_target_samples <- c(samples_ctrl, samples_dm1, samples_saha, samples_aso, samples_combo)

# ------------------------------------------------------------------------------
# 3. BUILD MASTER COUNT MATRIX
# ------------------------------------------------------------------------------
all_files <- list.files(path = FASTQ_DIR, pattern = "ReadsPerGene.out.tab$", full.names = TRUE)

# Map requested samples to actual file paths
my_files_all <- purrr::map_chr(all_target_samples, function(id) {
  match <- all_files[stringr::str_detect(all_files, id)]
  if(length(match) == 0) stop(paste("Missing file for:", id))
  return(match[1]) 
})

# Load and merge Column 2 (UNSTRANDED counts) for all 15 samples
count_matrix_all <- my_files_all %>%
  purrr::map(function(file) {
    read.table(file, skip = 4) %>%
      dplyr::select(V1, V2)
  }) %>%
  purrr::reduce(inner_join, by = "V1")

# Format rownames and colnames
rownames(count_matrix_all) <- count_matrix_all$V1
count_matrix_all <- count_matrix_all %>% dplyr::select(-V1)
colnames(count_matrix_all) <- all_target_samples

cat("Master Count Matrix Built. Dimensions:", dim(count_matrix_all), "\n")

# ------------------------------------------------------------------------------
# 4. DESeq2: PAIRWISE COMPARISONS
# ------------------------------------------------------------------------------
# Helper function to streamline pairwise DESeq2 execution
run_pairwise_deseq <- function(mat, sample_names, cond_levels) {
  subset_mat <- mat[, sample_names]
  condition <- factor(rep(cond_levels, each = 3), levels = cond_levels)
  colData <- data.frame(row.names = colnames(subset_mat), condition = condition)
  
  dds <- DESeqDataSetFromMatrix(subset_mat, colData, ~ condition)
  dds <- DESeq(dds)
  # Contrast maps to: c("condition", "Numerator (Treatment/Disease)", "Denominator (Control/Baseline)")
  res <- results(dds, contrast = c("condition", cond_levels[2], cond_levels[1]))
  return(res)
}

# A. DISEASE BASELINE: Healthy vs Affected (Positive LFC = Upregulated in DM1)
res_con_dm1 <- run_pairwise_deseq(count_matrix_all, c(samples_ctrl, samples_dm1), c("Healthy", "Affected"))

# B. ASO EFFICACY: Affected vs ASO (Positive LFC = Upregulated by Drug)
res_dm1_aso <- run_pairwise_deseq(count_matrix_all, c(samples_dm1, samples_aso), c("Affected", "ASO"))

# C. SAHA EFFICACY: Affected vs SAHA 
res_dm1_saha <- run_pairwise_deseq(count_matrix_all, c(samples_dm1, samples_saha), c("Affected", "SAHA"))

# D. COMBO EFFICACY: Affected vs Combo 
res_dm1_combo <- run_pairwise_deseq(count_matrix_all, c(samples_dm1, samples_combo), c("Affected", "Combo"))

print("Successfully generated 4 independent DESeq2 models.")

# ------------------------------------------------------------------------------
# 5. ISOLATE SIGNIFICANT DISEASE SIGNATURE
# ------------------------------------------------------------------------------
# Order and annotate baseline results
res_ordered <- res_con_dm1[order(res_con_dm1$padj), ]
res_df <- as.data.frame(res_ordered)
res_df$Symbol <- mapIds(org.Hs.eg.db, keys = rownames(res_df), column = "SYMBOL", 
                        keytype = "ENSEMBL", multiVals = "first")

write.csv(res_df, file.path(OUT_DIR, "Complete_Baseline_DGE(DMSO).csv"))

# Define the strict disease universe (padj < 0.05, |LFC| >= 1.5)
sig_degs <- subset(res_df, padj < 0.05 & abs(log2FoldChange) >= 1.5)
cat("Significant Disease DEGs found:", nrow(sig_degs), "\n")
write.csv(sig_degs, file.path(OUT_DIR, "Significant_Disease_DEGs(DMSO).csv"))

# ------------------------------------------------------------------------------
# 6. THERAPEUTIC RESCUE LOGIC (The Complete Universe)
# ------------------------------------------------------------------------------
# Consolidate LFCs across all comparisons for the target disease universe
rescue_lfc_table <- data.frame(
  gene_id      = rownames(res_con_dm1),
  LFC_Disease  = res_con_dm1$log2FoldChange,
  LFC_ASO      = res_dm1_aso$log2FoldChange,
  LFC_SAHA     = res_dm1_saha$log2FoldChange,
  LFC_Combo    = res_dm1_combo$log2FoldChange
) %>%
  filter(gene_id %in% rownames(sig_degs)) %>%
  left_join(sig_degs %>% rownames_to_column("gene_id") %>% dplyr::select(gene_id, Symbol), by = "gene_id")

# Calculate Rescue Percentages and Categorize (3-Bin Rule)
final_complete_dge <- rescue_lfc_table %>%
  mutate(
    Pct_ASO   = (LFC_ASO / -LFC_Disease) * 100,
    Pct_SAHA  = (LFC_SAHA / -LFC_Disease) * 100,
    Pct_Combo = (LFC_Combo / -LFC_Disease) * 100
  ) %>%
  mutate(across(starts_with("Pct_"), ~ case_when(
    is.na(.) ~ "NA",
    . >= 10  ~ "Rescue",
    . <= -10 ~ "Mis-rescue",
    TRUE     ~ "No change"
  ), .names = "Cat_{str_remove(.col, 'Pct_')}"))


print("--- Corrected Biological Universe ---")
print(table(final_complete_dge$Cat_Combo))

# ------------------------------------------------------------------------------
# 7. VISUALIZE RESCUE EFFICACY
# ------------------------------------------------------------------------------
# A. Clustered Bar Plot for Complete Distribution
plot_data <- final_complete_dge %>%
  pivot_longer(cols = starts_with("Cat_"), names_to = "Treatment", values_to = "Effect") %>%
  mutate(
    Treatment = factor(str_remove(Treatment, "Cat_"), levels = c("ASO", "SAHA", "Combo")),
    Effect = factor(Effect, levels = c("Rescue", "No change", "Mis-rescue", "NA"))
  )

p_bar <- ggplot(plot_data, aes(x = Treatment, fill = Effect)) +
  geom_bar(position = "dodge", color = "black", width = 0.7) +
  theme_minimal() +
  scale_fill_manual(values = c("Rescue" = "#2ca25f", "No change" = "#999999", 
                               "Mis-rescue" = "#de2d26", "NA" = "#fdbb84")) +
  labs(title = "Treatment Efficacy Across Disease Universe",
       subtitle = paste("Directional Analysis:", nrow(sig_degs), "Significant Genes"),
       y = "Gene Count", x = NULL) +
  theme(text = element_text(size = 14, face = "bold"), panel.grid.major.x = element_blank()) +
  geom_text(stat='count', aes(label=after_stat(count)), position=position_dodge(width=0.7), vjust=-0.5, size=3.5)

ggsave(file.path(OUT_DIR, "Final_Publication_Clustered_Bar.png"), p_bar, width = 8, height = 6, dpi = 300)

# B. Venn Diagram of Unique Gene Rescues
rescue_sets <- list(
  ASO   = final_complete_dge$gene_id[final_complete_dge$Pct_ASO >= 10],
  SAHA  = final_complete_dge$gene_id[final_complete_dge$Pct_SAHA >= 10],
  Combo = final_complete_dge$gene_id[final_complete_dge$Pct_Combo >= 10]
)

p_venn <- ggvenn(rescue_sets, fill_color = c("#377eb8", "#4daf4a", "#984ea3"), 
                 stroke_size = 0.5, set_name_size = 5, text_size = 4) +
  theme_minimal() +
  labs(title = "Uniqueness of Gene Rescue", subtitle = paste("Total Disease Signature:", nrow(sig_degs), "genes")) +
  theme(panel.grid = element_blank(), axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(OUT_DIR, "Rescued_Genes_Venn_Diagram.png"), p_venn, width = 8, height = 8, dpi = 300)

# ------------------------------------------------------------------------------
# 8. GENE ONTOLOGY (GO) ENRICHMENT: COMBO RESCUE
# ------------------------------------------------------------------------------
combo_rescued_ids <- final_complete_dge$gene_id[final_complete_dge$Pct_Combo >= 10]

go_combo <- enrichGO(
  gene          = combo_rescued_ids,
  OrgDb         = org.Hs.eg.db,
  keyType       = 'ENSEMBL',
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

p_go <- dotplot(go_combo, showCategory = 15) + 
  theme_minimal() +
  labs(title = "Biological Pathways Rescued by Combo", subtitle = "Genes reaching >= 10% rescue threshold") +
  theme(panel.grid = element_blank(), axis.line = element_line(color = "black"), text = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "Combo_Rescue_GO_Dotplot.png"), p_go, width = 10, height = 8, dpi = 300)
write.csv(as.data.frame(go_combo), file.path(OUT_DIR, "Combo_Rescue_GO_Results.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 9. OFF-TARGET TOXICITY ANALYSIS (The Stable Universe Method)
# ------------------------------------------------------------------------------
# A. Define the Stable Background (Genes lacking statistical certainty OR biological meaning in disease)
raw_baseline <- as.data.frame(res_con_dm1) %>% rownames_to_column("GeneID")
stable_baseline_degs <- raw_baseline %>%
  filter(is.na(padj) | padj >= 0.05 | abs(log2FoldChange) < 1.5)
stable_baseline_ids <- stable_baseline_degs$GeneID

# B. Extract highly significant drug-altered DEGs
prep_treatment_degs <- function(res_obj) {
  as.data.frame(res_obj) %>% rownames_to_column("GeneID") %>%
    filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1.5)
}

aso_degs   <- prep_treatment_degs(res_dm1_aso)
saha_degs  <- prep_treatment_degs(res_dm1_saha)
combo_degs <- prep_treatment_degs(res_dm1_combo)

# C. Isolate True Off-Target DEGs (Significant in treatment, but originated from stable background)
genes_off_target_aso   <- unique((aso_degs %>% filter(GeneID %in% stable_baseline_ids))$GeneID)
genes_off_target_saha  <- unique((saha_degs %>% filter(GeneID %in% stable_baseline_ids))$GeneID)
genes_off_target_combo <- unique((combo_degs %>% filter(GeneID %in% stable_baseline_ids))$GeneID)

# D. Visualise Shared Toxicity (UpSet Plot)
off_target_list <- list(ASO = genes_off_target_aso, SAHA = genes_off_target_saha, Combo = genes_off_target_combo)

png(file.path(OUT_DIR, "Off_Target_DEGs_UpSet1.png"), width = 10, height = 8, units = "in", res = 300)
upset(fromList(off_target_list), nsets = 3, nintersects = 7, main.bar.color = "black",
      sets.bar.color = c("#377EB8", "#4DAF4A", "#984EA3"), text.scale = c(1.5, 1.5, 1.2, 1.2, 1.5, 1.3),
      order.by = "freq", mainbar.y.label = "Number of Off-Target DEGs", sets.x.label = "Total Off-Target DEGs per Treatment")
grid.text("Off-Target Effect on DEG Analysis", x = 0.65, y = 0.95, gp = gpar(fontsize = 15, fontface = "bold"))
dev.off() 

# E. Extract and Export the Toxicity Intersections
extracted_degs <- data.frame(
  Category = c(
    rep("All_Three_Shared", length(Reduce(intersect, off_target_list))),
    rep("Combo_Unique", length(setdiff(genes_off_target_combo, union(genes_off_target_saha, genes_off_target_aso)))),
    rep("ASO_and_Combo", length(setdiff(intersect(genes_off_target_aso, genes_off_target_combo), genes_off_target_saha))),
    rep("SAHA_and_Combo", length(setdiff(intersect(genes_off_target_saha, genes_off_target_combo), genes_off_target_aso)))
  ),
  GeneID = c(
    Reduce(intersect, off_target_list),
    setdiff(genes_off_target_combo, union(genes_off_target_saha, genes_off_target_aso)),
    setdiff(intersect(genes_off_target_aso, genes_off_target_combo), genes_off_target_saha),
    setdiff(intersect(genes_off_target_saha, genes_off_target_combo), genes_off_target_aso)
  )
)

write.csv(extracted_degs, file.path(OUT_DIR, "Extracted_OffTarget_DEGs_Intersections.csv"), row.names = FALSE)
cat("Pipeline Complete. Extracted Off-Target Gene IDs saved.\n")