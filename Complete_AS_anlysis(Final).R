# ==============================================================================
# SCRIPT: Myotonic Dystrophy (DM1) Alternative Splicing & Rescue Pipeline
# DESCRIPTION: Comprehensive analysis of rMATS splicing data, including baseline 
#              disease characterization, visualization of key biomarkers, 
#              therapeutic rescue efficacy (ASO, SAHA, Combo), and off-target 
#              toxicity evaluation.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. SETUP & LIBRARIES
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  # Core Data Manipulation & Plotting
  library(tidyverse)    # Covers dplyr, tidyr, ggplot2, stringr, etc.
  library(patchwork)
  library(scales)
  
  # Genomics & Splicing
  library(maser)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(GenomicFeatures)
  library(GenomicRanges)
  
  # Visualization Extras
  library(ggvenn)
  library(UpSetR)
  
  # Pathway Analysis
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
})

# ------------------------------------------------------------------------------
# 2. CONFIGURATION & FILE PATHS
# ------------------------------------------------------------------------------
# Define base directories here so the script is easily reproducible on other machines
DATA_DIR <- "C:/Users/USER/Documents/Berglund/Phase_two/Phase 2/data/AS-data"
PLOT_DIR <- "C:/Users/USER/Documents/Berglund/For_Presentation/Phase_1"
MASER_DIR <- "Z:\\common\\Tilesanmi_SFD\\First_Project\\Phase_1\\Taye_STAR\\BAM_Files\\AS-Resul"

# Define specific rMATS file paths
file_baseline <- file.path(DATA_DIR, "SE.MATS.JCEC.txt_DM05_vs_Control.txt")
file_aso      <- file.path(DATA_DIR, "SE.MATS.JCEC_DM05_MB_DMSO_vs_DM05_MB_ASO_15.txt")
file_saha     <- file.path(DATA_DIR, "SE.MATS.JCEC.txt_DM05_vs_SAHA.txt")
file_combo    <- file.path(DATA_DIR, "SE.MATS.JCEC.txt_DM05_VS_AS0+SAHA.txt") 

# Set working directory to data location for maser initialization

# ------------------------------------------------------------------------------
# 3. GLOBAL SPLICING ALTERATIONS (DM1 vs CONTROL)
# ------------------------------------------------------------------------------
# Initialize Maser object
Disease <- maser(MASER_DIR, c("DM", "Control"), ftype = "JCEC")

# Extract summary tables for all event types
df_se   <- summary(Disease, type = "SE")
df_RI   <- summary(Disease, type = "RI")
df_MXE  <- summary(Disease, type = "MXE")
df_A3SS <- summary(Disease, type = "A3SS")
df_A5SS <- summary(Disease, type = "A5SS")

# Helper function to filter for significant events (FDR <= 0.05, |dPSI| >= 10%)
filter_sig <- function(df) {
  df %>% filter(FDR <= 0.05 & abs(IncLevelDifference) >= 0.1)
}

# Apply strict filtering
sig_se   <- filter_sig(df_se)
sig_RI   <- filter_sig(df_RI)
sig_MXE  <- filter_sig(df_MXE)
sig_A3SS <- filter_sig(df_A3SS)
sig_A5SS <- filter_sig(df_A5SS)

# Prepare data for the Global Summary Bar Plot
plot_data_summary <- data.frame(
  Category = c("Skipped Exon (SE)", "Mutually Exclusive (MXE)", 
               "Retained Intron (RI)", "Alt 5' SS", "Alt 3' SS"),
  Count = c(nrow(sig_se), nrow(sig_MXE), nrow(sig_RI), nrow(sig_A5SS), nrow(sig_A3SS))
)

# Plot: Global Summary
final_summary_plot <- ggplot(plot_data_summary, aes(x = reorder(Category, -Count), y = Count, fill = Category)) +
  geom_bar(stat = "identity", color = "black", linewidth = 1, width = 0.7) +
  geom_text(aes(label = Count), vjust = -0.5, fontface = "bold", size = 5) + 
  scale_fill_manual(values = c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00")) +
  theme_classic(base_size = 16) +
  labs(title = "Global Splicing Alterations: DM1 vs Control",
       subtitle = "Significant Events (FDR < 0.05, |ΔPSI| > 10%)",
       x = "Splicing Category", y = "Number of Significant Events") +
  theme(legend.position = "none",
        axis.text.x = element_text(face = "bold", angle = 15, hjust = 1),
        plot.title = element_text(face = "bold", size = 20))

ggsave(file.path(PLOT_DIR, "Correct_DM1_Global_Summary.png"), final_summary_plot, width = 10, height = 7, dpi = 300)

# Plot: PCA of Splicing Events
p_pca <- pca(Disease, type = "SE") + 
  theme_classic(base_size = 18) + 
  theme(
    axis.line = element_line(linewidth = 1.2, color = "black"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black", face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  geom_point(size = 6)

ggsave(file.path(PLOT_DIR, "PCA_Splicing_DM1.png"), plot = p_pca, width = 8, height = 6, dpi = 300, bg = "white")

# ------------------------------------------------------------------------------
# 4. TARGETED BIOMARKER VISUALIZATION (INSR, BIN1, MBNL1, SYNE1)
# ------------------------------------------------------------------------------
# Annotate specific exons using UCSC hg38 database
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
all_exons <- exons(txdb, columns = c("exon_id", "gene_id"))

# Convert rMATS coordinates (0-based) to UCSC ranges (1-based)
se_coords <- sig_se %>% separate(exon_target, into = c("Start", "End"), sep = "-", convert = TRUE)
query_ranges <- GRanges(
  seqnames = as.character(sig_se$Chr),
  ranges = IRanges(start = se_coords$Start + 1, end = se_coords$End),
  strand = as.character(sig_se$Strand)
)

# Map database IDs to significant SE events
hits <- findOverlaps(query_ranges, all_exons)
sig_se$Database_Exon_ID <- NA
sig_se$Database_Exon_ID[queryHits(hits)] <- all_exons$exon_id[subjectHits(hits)]

# Define target biomarkers and their UCSC IDs
target_exons_manual <- list(
  "INSR"  = list(id = "794738", label = "INSR: Exon 11"),
  "BIN1"  = list(id = "136027", label = "BIN1: Exon 11"),
  "MBNL1" = list(id = "171836", label = "MBNL1: Exon 5"),
  "SYNE1" = list(id = "330615", label = "SYNE1: Exon 137")
)

# Plotting function for individual biomarker panels
make_manual_plot <- function(gene_name) {
  target_id    <- target_exons_manual[[gene_name]]$id
  target_label <- target_exons_manual[[gene_name]]$label
  
  event_data <- sig_se %>% filter(geneSymbol == gene_name, Database_Exon_ID == target_id)
  if (nrow(event_data) == 0) return(NULL)
  
  parse_psi <- function(x) as.numeric(unlist(strsplit(as.character(x), ","))) %>% na.omit()
  psi_dm   <- parse_psi(event_data$PSI_1[1])
  psi_ctrl <- parse_psi(event_data$PSI_2[1])
  
  plot_df <- data.frame(
    Group = factor(c(rep("DM", length(psi_dm)), rep("Control", length(psi_ctrl))), 
                   levels = c("Control", "DM")),
    PSI = c(psi_dm, psi_ctrl)
  )
  sum_df <- plot_df %>% group_by(Group) %>% summarise(m = mean(PSI), s = sd(PSI))
  
  ggplot(sum_df, aes(x = Group, y = m, fill = Group)) +
    geom_bar(stat = "identity", color = "black", linewidth = 1.2, width = 0.6) +
    geom_errorbar(aes(ymin = m - s, ymax = m + s), width = 0.2, linewidth = 1) +
    geom_jitter(data = plot_df, aes(x = Group, y = PSI), width = 0.1, size = 4, color = "black") +
    scale_fill_manual(values = c("DM" = "#E41A1C", "Control" = "#377EB8")) +
    theme_classic(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(face = "bold", color = "black"),
      legend.position = "none",
      axis.line = element_line(linewidth = 1)
    ) +
    labs(title = target_label, y = "Inclusion Level (PSI)", x = "") +
    scale_y_continuous(labels = percent, limits = c(0, 1), expand = c(0,0))
}

# Assemble Patchwork Grid for Biomarkers
biomarker_grid <- (make_manual_plot("INSR") | make_manual_plot("BIN1")) / 
  (make_manual_plot("MBNL1") | make_manual_plot("SYNE1")) + 
  plot_annotation(
    title = "Key Splicing Biomarkers in Myotonic Dystrophy (DM1)",
    theme = theme(plot.title = element_text(size = 22, face = "bold", hjust = 0.5))
  )

ggsave(file.path(PLOT_DIR, "Key_missplicing_events_SE.png"), biomarker_grid, width = 12, height = 10, dpi = 300)

# ------------------------------------------------------------------------------
# 5. THERAPEUTIC RESCUE & MISRESCUE ANALYSIS
# ------------------------------------------------------------------------------
# Helper to average comma-separated PSI strings robustly
get_mean <- function(x) {
  sapply(strsplit(as.character(x), ","), function(vals) {
    num_vals <- as.numeric(vals)
    if(all(is.na(num_vals))) return(NA)
    mean(num_vals, na.rm = TRUE)
  })
}

# Helper to create an immutable coordinate fingerprint for tracking events across runs
make_uni_id <- function(df) {
  df %>% mutate(New_Uni_ID = paste0(GeneID, chr, strand, exonStart_0base, exonEnd, 
                                    upstreamES, upstreamEE, downstreamES, downstreamEE))
}

# Load raw files with fingerprints applied
d_con <- read.delim(file_baseline) %>% make_uni_id()
d_aso <- read.delim(file_aso) %>% make_uni_id()
d_sah <- read.delim(file_saha) %>% make_uni_id()
d_com <- read.delim(file_combo) %>% make_uni_id()

# Isolate the Disease Universe (Significant baseline events)
universe <- d_con %>%
  filter(FDR < 0.05 & abs(IncLevelDifference) >= 0.1) %>%
  mutate(PSI_DM1 = get_mean(IncLevel1), PSI_Healthy = get_mean(IncLevel2)) %>%
  dplyr::select(New_Uni_ID, GeneID, geneSymbol, PSI_DM1, PSI_Healthy)

# Map treatments to the universe and calculate rescue metrics
mapped_data <- universe %>%
  left_join(d_aso %>% dplyr::select(New_Uni_ID, IncLevel2_ASO = IncLevel2), by = "New_Uni_ID") %>%
  mutate(PSI_ASO = get_mean(IncLevel2_ASO)) %>%
  left_join(d_sah %>% dplyr::select(New_Uni_ID, IncLevel2_SAHA = IncLevel2), by = "New_Uni_ID") %>%
  mutate(PSI_SAHA = get_mean(IncLevel2_SAHA)) %>%
  left_join(d_com %>% dplyr::select(New_Uni_ID, IncLevel2_Combo = IncLevel2), by = "New_Uni_ID") %>%
  mutate(PSI_Combo = get_mean(IncLevel2_Combo)) %>%
  mutate(
    Pct_ASO   = (1 - ((PSI_Healthy - PSI_ASO) / (PSI_Healthy - PSI_DM1))) * 100,
    Pct_SAHA  = (1 - ((PSI_Healthy - PSI_SAHA) / (PSI_Healthy - PSI_DM1))) * 100,
    Pct_Combo = (1 - ((PSI_Healthy - PSI_Combo) / (PSI_Healthy - PSI_DM1))) * 100
  )

# Categorize using the 3-Bin Rule (10% threshold)
final_data <- mapped_data %>%
  mutate(across(starts_with("Pct_"), ~ case_when(
    is.na(.) ~ "Not Detected in Treatment",
    . >= 10  ~ "Rescue",
    . <= -10 ~ "Misrescue",
    TRUE     ~ "No Change"
  ), .names = "Cat_{str_remove(.col, 'Pct_')}"))

# Plot: Splicing Rescue Distribution (Clustered Bar)
plot_data_rescue <- final_data %>%
  pivot_longer(cols = starts_with("Cat_"), names_to = "Treatment", values_to = "Effect") %>%
  mutate(Treatment = factor(str_remove(Treatment, "Cat_"), levels = c("ASO", "SAHA", "Combo")),
         Effect = factor(Effect, levels = c("Rescue", "No Change", "Misrescue", "Not Detected in Treatment")))

p_distribution <- ggplot(plot_data_rescue, aes(x = Treatment, fill = Effect)) +
  geom_bar(position = "dodge", color = "black", width = 0.7) +
  theme_minimal() +
  scale_fill_manual(
    values = c("Rescue" = "#2ca25f", "No Change" = "#999999", 
               "Misrescue" = "#de2d26", "Not Detected in Treatment" = "#fdbb84"),
    drop = FALSE 
  ) +
  labs(title = "Splicing Rescue Distribution",
       subtitle = paste("Directional Analysis:", nrow(universe), "Significant Events"),
       y = "Gene Count", x = NULL) +
  theme(text = element_text(size = 14, face = "bold"), panel.grid.major.x = element_blank()) +
  geom_text(stat='count', aes(label=after_stat(count)), position=position_dodge(width=0.7), vjust=-0.5, size=3.5)

ggsave(file.path(PLOT_DIR, "Final_Publication_Clustered_Bar_Splicing.png"), p_distribution, width = 8, height = 6, dpi = 300)

# Plot: Overlap of Rescued Events (Venn Diagram)
venn_list <- list(
  ASO   = final_data %>% filter(Cat_ASO == "Rescue") %>% pull(New_Uni_ID) %>% unique(),
  SAHA  = final_data %>% filter(Cat_SAHA == "Rescue") %>% pull(New_Uni_ID) %>% unique(),
  Combo = final_data %>% filter(Cat_Combo == "Rescue") %>% pull(New_Uni_ID) %>% unique()
)

venn_plot <- ggvenn(venn_list, show_elements = FALSE, fill_color = c("#3182bd", "#2ca25f", "#756bb1"), 
                    fill_alpha = 0.6, stroke_size = 0.5, set_name_size = 5) +
  labs(title = "Overlap of Rescued Splicing Events", subtitle = "Events Classified as 'Rescue' (>= 10% Recovery)") +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5), plot.subtitle = element_text(hjust = 0.5, size = 12))

ggsave(file.path(PLOT_DIR, "Venn_Rescued_Events.png"), plot = venn_plot, width = 6, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# 6. GENE ONTOLOGY (GO) ENRICHMENT: HIGH EFFICACY RESCUE
# ------------------------------------------------------------------------------
# Identify highly rescued genes in the Combo treatment (>50% rescue)
baseline <- d_con %>%
  filter(FDR < 0.05) %>%
  mutate(PSI_Con = get_mean(IncLevel1), PSI_DM1 = get_mean(IncLevel2), Delta_Disease = PSI_DM1 - PSI_Con) %>%
  dplyr::select(New_Uni_ID, geneSymbol, PSI_Con, PSI_DM1, Delta_Disease)

combo_rescued_genes <- d_com %>%
  filter(FDR < 0.05) %>%
  mutate(PSI_Trt = get_mean(IncLevel2)) %>%
  dplyr::select(New_Uni_ID, geneSymbol, PSI_Trt) %>%
  inner_join(baseline, by = c("New_Uni_ID", "geneSymbol")) %>%
  mutate(
    Delta_Rescue = PSI_Trt - PSI_DM1,
    is_rescued = sign(Delta_Rescue) != sign(Delta_Disease) & (abs(Delta_Rescue) >= 0.5 * abs(Delta_Disease))
  ) %>%
  filter(is_rescued == TRUE) %>% pull(geneSymbol) %>% unique()

# Run Enrichment if genes are present
if(length(combo_rescued_genes) > 0) {
  go_results_combo <- enrichGO(gene = combo_rescued_genes, OrgDb = org.Hs.eg.db, 
                               keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.1)
  
  go_dotplot <- dotplot(go_results_combo, showCategory = 20) + 
    ggtitle("Biological Processes Enriched in Combo-Rescued Events") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14))
  
  ggsave(file.path(PLOT_DIR, "Combo_Rescue_GO_Dotplot.png"), plot = go_dotplot, width = 8, height = 8, dpi = 300)
} else {
  print("No genes found for GO Enrichment.")
}

# ------------------------------------------------------------------------------
# 7. OFF-TARGET TOXICITY ANALYSIS
# ------------------------------------------------------------------------------
# Helper function to strictly isolate events significantly altered by drug treatment
process_rmats_toxicity <- function(filepath) {
  read.table(filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
    mutate(Composite_ID = paste(chr, strand, exonStart_0base, exonEnd, sep = "_")) %>%
    filter(FDR < 0.05 & abs(IncLevelDifference) >= 0.1)
}

# Define the "Stable Background Universe" using the UNFILTERED baseline data.
# Inverse logic captures events lacking statistical certainty OR biological meaning.
stable_baseline_ids <- read.table(file_baseline, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  mutate(Composite_ID = paste(chr, strand, exonStart_0base, exonEnd, sep = "_")) %>%
  filter(FDR >= 0.05 | abs(IncLevelDifference) < 0.1) %>% 
  pull(Composite_ID) %>% 
  unique()

# Extract True Off-Target Toxicity (Significant in drug, but originated from the stable background)
events_off_target_aso   <- process_rmats_toxicity(file_aso) %>% filter(Composite_ID %in% stable_baseline_ids) %>% pull(Composite_ID) %>% unique()
events_off_target_saha  <- process_rmats_toxicity(file_saha) %>% filter(Composite_ID %in% stable_baseline_ids) %>% pull(Composite_ID) %>% unique()
events_off_target_combo <- process_rmats_toxicity(file_combo) %>% filter(Composite_ID %in% stable_baseline_ids) %>% pull(Composite_ID) %>% unique()

# Plot: Shared Toxicity (UpSet Plot)
off_target_list <- list(ASO = events_off_target_aso, SAHA = events_off_target_saha, Combo = events_off_target_combo)

png(file.path(PLOT_DIR, "Off_Target_Toxicity_UpSet.png"), width = 800, height = 600, res = 120)
upset(fromList(off_target_list), order.by = "freq", 
      mainbar.y.label = "Shared Off-Target Events", sets.x.label = "Total Toxicity per Treatment",
      text.scale = c(1.3, 1.3, 1, 1, 1.5, 1))
dev.off()