### Transcriptomic Analysis of Therapeutic Rescue in Myotonic Dystrophy Type 1 (DM1)
This repository contains bioinformatics pipelines developed during a virtual research training program at the Berglund Lab, University at Albany. The project evaluates the therapeutic efficacy and transcriptomic impact of Antisense Oligonucleotides (ASO), Suberoylanilide Hydroxamic Acid (SAHA), and combination treatments on patient-derived Myotonic Dystrophy type 1 (DM1) myoblasts.
##Project Overview
Myotonic Dystrophy type 1 is characterized by widespread alternative splicing (AS) defects and downstream differential gene expression (DGE). This pipeline processes raw RNA-seq outputs to quantify how well different therapeutic interventions rescue the disease signature, while strictly monitoring for drug-induced off-target toxicity.


## Repository Contents
#Script 
Complete end to end differential gene expression analysis and alternative splicing analysis scripts 
#ppt
A ppt summarizing all my findings from the beginning to end of the project 


## Methodology
# 1. Robust Therapeutic Rescue Logic 
To quantify therapeutic efficacy, a strictly defined "Disease Universe" is established for both AS and DGE data. Treatment effects are mapped to this universe and binned into three distinct functional categories:

Rescue: 10% or greater reversal toward the healthy phenotype.

Mis-rescue: 10% or greater exacerbation of the disease phenotype.

No Change: Negligible transcriptomic movement.

# 2. Off-Target Toxicity Profiling
Identifying true drug-induced toxicity requires bypassing statistical thresholding artifacts. This pipeline utilizes a top-down Stable Background Universe approach:

Define the Baseline: First, we isolate transcriptomic events (genes or splicing coordinates) that definitively did not change in the disease state (e.g., FDR >= 0.05 or sub-threshold effect sizes).

Track Drug Impact: Second, we filter for events within this stable universe that became significantly mis-regulated only after drug administration.

Biological Verification: This ensures that reported off-target effects represent genuine de novo drug toxicity rather than borderline disease-associated events pushed over an arbitrary statistical threshold.


## Tool Stack & Dependencies
The complete analytical pipeline spans upstream execution on the UAlbany HPC cluster and downstream statistical modeling in R.

# Upstream Processing (UAlbany HPC)

FastQC: Raw read quality control and visualization.

fastp: High-performance adapter trimming, quality filtering, and read pruning.

STAR (Spliced Transcripts Alignment to a Reference): Rapid, splice-aware alignment of RNA-seq reads to the human reference genome and generation of gene count matrices.

rMATS (replicate Multivariate Analysis of Transcript Splicing): Detection and quantification of differential alternative splicing events.

# Downstream Analysis (R Environment)
Core & Data Wrangling: tidyverse, stringr

Transcriptomics & Splicing: DESeq2, maser, TxDb.Hsapiens.UCSC.hg38.knownGene, GenomicFeatures, GenomicRanges

Visualization: ggplot2, patchwork, ggvenn, UpSetR, grid, scales

Ontology & Enrichment: clusterProfiler, org.Hs.eg.db, AnnotationDbi, enrichplot

## Acknowledgments
Developed during virtual research training at the University at Albany, supporting ongoing investigations in the Berglund Lab into targeted therapeutics for DM1. Special thanks to Porama Hoque and Claudia Lennon for their support. 
