library(TSRexploreR)
library(tidyverse)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)

# Read in sample sheet
samples <- read_tsv("samples_ctss.txt")

# Create tsrexplorer object
exp <- tsr_explorer(sample_sheet = samples, 
                    genome_assembly = BSgenome.Scerevisiae.UCSC.sacCer3,
                    genome_annotation = TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)

# Import CTSSs
exp <- tss_import(exp, sample_sheet = samples)

# Format counts
exp <- format_counts(exp, data_type = "tss")

# Annotate TSSs
exp <- annotate_features(exp, data_type = "tss", feature_type = "transcript", 
                         upstream = 250, downstream = 100)

# Explore thresholds (Fig. 1B)
threshold_data <- explore_thresholds(exp, steps = 1, max_threshold = 50, 
                                     samples = "all", use_normalized = FALSE)

plot_threshold_exploration(threshold_data, ncol = 4, point_size = .5) +
  scale_color_viridis_c() 

# Genomic distribution (Fig. 1C)
tss_distribution <- genomic_distribution(exp, data_type = "tss", threshold = 10)

plot_genomic_distribution(tss_distribution) +
  scale_fill_viridis_d(direction = -1, name = "Annotation")

# Detected features (Fig. 1D)
features <- detect_features(exp, data_type = "tss", threshold = 10)

plot_detected_features(features) +
  scale_fill_viridis_d(direction = -1)

# Normalize counts
exp <- normalize_counts(exp, data_type = "tss", method = "deseq2", 
                        threshold = 10, n_samples = 1)

# Correlation (Fig. 1E)
plot_correlation(
  exp, data_type = "tss", 
  font_size = 12,
  use_normalized = TRUE, 
  cluster_samples = TRUE, 
  correlation_metric = "pearson",
  heatmap_colors = viridis::viridis(256)
)

# Tracks (Fig. 1F)
gene_tracks(exp, feature_name = "YLR081W", promoter_only = FALSE, 
            samples = c(TSS = "YPD_1", TSS = "galactose_1"), 
            ymax = 10000, tss_colors = viridis::viridis(4), 
            use_normalized = TRUE, axis_scale = 1)

# Density plot (Fig. 1G)
plot_density(exp, data_type = "tss", samples = "YPD_1", threshold = 10)

# Heatmap (Fig. 1H)
tss_matrix <- tss_heatmap_matrix(exp, use_normalized = TRUE, samples = "YPD_1")

plot_heatmap(tss_matrix, background_color = "white")

# Max UTR length (Fig. 1I)
max_utr_lengths <- max_utr(exp, samples = "YPD_1", threshold = 10, 
                           feature_type = "transcriptId")

plot_max_utr(max_utr_lengths)

####################
### TSR analysis ###
####################

# Cluster TSSs
exp <- tss_clustering(exp, max_distance = 25, max_width = 100, threshold = 10)

# Associate TSSs with TSRs
exp <- associate_with_tsr(exp)

exp <- format_counts(exp, data_type = "tsr")

# TSR metrics (Fig. 2A)
exp <- tsr_metrics(exp)

plot_tsr_metric(exp, tsr_metrics = c("shape_index", "iqr_width", "peak_balance"), 
                log2_transform = FALSE, ncol = 3, 
                samples = c("YPD_1", "arrest_1", "diauxic_shift_1", "dna_damage_1",
                            "galactose_1", "glucose_1", "H2O2_1", "heat_shock_1",
                            "NaCl_1"))

# Sequence analysis
exp <- mark_dominant(exp, data_type = "tss", threshold = 10, use_normalized = FALSE)

# Sequence logos (Fig. 2B-C)
conditions <- list(quantile_by = "iqr_width",
                   n_quantiles = 5,
                   quantile_direction = "descending",
                   filter = "score > 25")

seqs <- tss_sequences(exp, dominant = TRUE, 
                      data_conditions = conditions, samples = "YPD_1")

plot_sequence_logo(seqs, ncol = 2)

# Color map (Fig. 2D)
plot_sequence_colormap(seqs, ncol = 4)

# Dinucleotide frequencies (Fig. 2E)
frequencies <- dinucleotide_frequencies(exp, dominant = TRUE, samples = "all")

plot_dinucleotide_frequencies(frequencies, ncol = 4) +
  scale_fill_viridis_c()

# Annotate TSRs
exp <- annotate_features(exp, data_type = "tsr", feature_type = "transcript",
                         upstream = 250, downstream = 100)

# Correlation
plot_correlation(
  exp, data_type = "tsr", 
  font_size = 12,
  use_normalized = TRUE,
  cluster_samples = TRUE, 
  correlation_metric = "pearson",
  heatmap_colors = viridis::viridis(100)
)

#################################
### Differential TSR analysis ###
#################################

# Perform DE analysis
exp <- fit_de_model(exp, data_type = "tsr", formula = ~condition, method = "deseq2")

treatments <- filter(exp@meta_data$sample_sheet, condition != "YPD") %>%
  dplyr::select(condition) %>%
  distinct() %>%
  unlist()

# exp@diff_features$TSRs$results <- NULL

# Loop to compare all treatments to YPD
for(i in treatments) { 
  exp <- differential_expression(
    exp, data_type = "tsr", 
    comparison_name = str_c(i, "_vs_YPD"),
    comparison_type = "contrast",
    comparison = c("condition", i, "YPD"))
  }

# Plot number of DE TSRs (Fig. 2F)
plot_num_de(exp, data_type = "tsr", de_comparisons = "all",
            log2fc_cutoff = 1, fdr_cutoff = 0.05, 
            keep_unchanged = FALSE) +
  scale_fill_viridis_d() +
  theme_bw()

# MA plot (Fig. 2G)
plot_ma(exp, data_type = "tsr", de_comparisons = "heat_shock_vs_YPD", size = 0.5) +
  scale_color_viridis_d() +
  theme_bw()

# Volcano plot (Fig. 2H)
plot_volcano(exp, data_type = "tsr", de_comparisons = "heat_shock_vs_YPD", size = 0.5) +
  scale_color_viridis_d() +
  theme_bw()

# Annotate DE TSRs
exp <- annotate_features(exp, data_type = "tsr_diff", feature_type = "transcript",
                         upstream = 250, downstream = 100)

# GO analysis (Fig. 2I)
enrichment_data <- export_for_enrichment(exp, data_type = "tsr", 
                                         de_comparisons = "heat_shock_vs_YPD",
                                         keep_unchanged = FALSE) %>%
  filter(simple_annotations == "Promoter")

library("clusterProfiler")
library("org.Sc.sgd.db")

go_enrichment <- compareCluster(
  geneId ~ sample + de_status,
  data = enrichment_data,
  fun = "enrichGO",
  OrgDb = "org.Sc.sgd.db",
  pAdjustMethod = "fdr",
  ont = "BP",
  keyType="ENSEMBL"
)

dotplot(go_enrichment, font.size = 12, showCategory = 25) +
  scale_color_viridis_c() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
