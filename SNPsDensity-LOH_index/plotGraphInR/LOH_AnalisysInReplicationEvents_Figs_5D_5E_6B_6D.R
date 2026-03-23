###################################################################################################
# R script for LOH Analysis (Genomic Instability)
#
# Overview:
# The Loss of Heterozygosity (LOH) index analysis was calculated in the SNP analysis 
# and variant calling pipeline from the SNP VCF files. This same script was used to 
# calculate the loss of heterozygosity in the region of origins and terminations 
# ("origins", "IZD", "terminations", "TZD").
#
# Data Input: 1 BedGraph containing LOH Index
# Methodology: 
#   - Robust Stats: Shapiro-Wilk/Levene -> ANOVA/Tukey OR Kruskal-Wallis/Dunn
#   - Flexible Alpha: Adjustable threshold for normality and significance.
#   - Visuals: Metaplot + Boxplot with geom_signif and CSV Export.
###################################################################################################

# =============================
# Step 1: 📦 Load required packages
# =============================
required_packages <- c("data.table", "ggplot2", "cowplot", "ggpubr", "car", "tools", "rstatix", "ggsignif")
bioc_packages <- c("GenomicRanges")

check_install_packages <- function(pkg_list, is_bioc = FALSE) {
  for (pkg in pkg_list) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (is_bioc) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
      } else {
        install.packages(pkg)
      }
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

check_install_packages(required_packages)
check_install_packages(bioc_packages, is_bioc = TRUE)

# =============================
# Step 2: 📁 Input variables & Setup
# =============================
SUBDIR <- "LOH_Analisys/"
base_dir <- file.path("~/Documents/DNAscent-Origins/dataAnalisys/Profile", SUBDIR)

if (!dir.exists(base_dir)) {
  dir.create(base_dir, recursive = TRUE)
}
setwd(base_dir)

# --- 💡 STATISTICAL SENSITIVITY ADJUSTMENT ---
# If you want to accept larger p-values (e.g., trends), change to 0.10
ALPHA_THRESHOLD <- 0.0503 

# --- INPUT FILES ---
bed_files <- c("control.bed", "terminations.bed", "miniTerminationZones.bed")  
bedgraph_file <- "average_LOH_consensus.bedgraph"

# Parameters
bin_size <- 20          
region_size <- 2000 
nbins <- region_size / bin_size
Experiment  <- "LOH Index (Median)" 

# =============================
# Step 3: 📥 Load bedGraph & Global Stats
# =============================
message("Loading LOH bedGraph file...")
bg_data <- fread(bedgraph_file, header = FALSE, select = 1:4)
setnames(bg_data, c("chr", "start", "end", "score"))

global_mean <- mean(bg_data$score, na.rm = TRUE)
bg.gr <- GRanges(seqnames = bg_data$chr, 
                 ranges = IRanges(start = bg_data$start + 1, end = bg_data$end), 
                 score = bg_data$score)

# =============================
# Step 4: 🛠️ Processing function
# =============================
process_bed_file <- function(bed_file, bg.gr, nbins, region_size) {
  bed <- fread(bed_file, header = FALSE)
  setnames(bed, 1:3, c("chr", "start", "end")) 
  
  gr <- GRanges(seqnames = bed$chr, ranges = IRanges(start = bed$start + 1, end = bed$end))
  gr.centered <- resize(gr, width = region_size, fix = "center")
  gr.bins <- unlist(tile(gr.centered, n = nbins))
  
  ov <- findOverlaps(gr.bins, bg.gr)
  hits_dt <- data.table(queryHits = queryHits(ov), subjectHits = subjectHits(ov))
  hits_dt[, bin_id := ((queryHits - 1) %% nbins) + 1]
  hits_dt[, region_id := ((queryHits - 1) %/% nbins) + 1]
  hits_dt[, val := bg.gr$score[subjectHits]]
  
  stats_per_bin <- hits_dt[, .(score = median(val, na.rm = TRUE)), by = .(region_id, bin_id)]
  all_bins_dt <- data.table(region_id = rep(1:length(gr), each = nbins), bin_id = rep(1:nbins, length(gr)))
  bins_dt <- merge(all_bins_dt, stats_per_bin, by = c("region_id", "bin_id"), all.x = TRUE)
  
  prof <- bins_dt[, .(mean_signal = median(score, na.rm = TRUE)), by = bin_id]
  file_type <- tools::file_path_sans_ext(basename(bed_file))
  prof$type <- file_type
  bins_dt$type <- file_type
  
  return(list(profile = prof, raw = bins_dt))
}

# =============================
# Step 5: 📊 Process Data
# =============================
all_results <- lapply(bed_files, process_bed_file, bg.gr = bg.gr, nbins = nbins, region_size = region_size)
df <- rbindlist(lapply(all_results, `[[`, "profile"))
df_raw <- rbindlist(lapply(all_results, `[[`, "raw"))

# =============================
# Step 6: 🏷️ Rename and Factor
# =============================
rename_samples <- function(type_name) {
  switch(type_name,
         "control" = "without \n terminations",
         "terminations" = "termination",
         "miniTerminationZones" = "TZD" ,
         type_name)
}
df[, type := sapply(type, rename_samples)]
df_raw[, type := sapply(type, rename_samples)]
levels_order <- c("without \n terminations", "termination", "TZD" )
df$type <- factor(df$type, levels = levels_order)
df_raw$type <- factor(df_raw$type, levels = levels_order)

# =============================
# Step 7: 🎨 Colors and Position
# =============================
color_group <- c(
  "without \n terminations" = rgb(44, 123, 182,  maxColorValue = 255), 
  "termination"              = rgb(100, 100, 255, maxColorValue = 255),
  "TZD"                      = rgb(255, 25, 0,    maxColorValue = 255)   
)

positions <- seq(-region_size/2 + bin_size/2, region_size/2 - bin_size/2, by = bin_size)
df$position <- positions[df$bin_id]

# Label for the reference line in the plot
genome_label <- "Genome\nMean"

# =============================
# Step 8: 📈 Plot A - Metaplot
# =============================
df_plot <- df[mean_signal > 0 & !is.na(mean_signal)]

pA <- ggplot(df_plot, aes(x = position, y = mean_signal, color = type)) +
  geom_line(linewidth = 1) + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size=0.2) +
  
  # Adding Genome Mean with line break in the legend
  geom_hline(aes(yintercept = global_mean, linetype = genome_label), 
             color = "red", size = 0.4, alpha = 0.7) +
  
  scale_linetype_manual(name = "", values = setNames("dashed", genome_label)) +
  scale_x_continuous(breaks = c(-region_size/2, 0, region_size/2), labels = c("S", "C", "E")) +
  scale_color_manual(values = color_group) +
  labs(x = NULL, y = "LOH index", color = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    panel.background = element_rect(color = "black", fill = NA), 
    legend.position = "none")
pA
# =============================
# Step 9: 🔍 Robust Statistical Analysis
# =============================
message("--- Statistical Analysis (Alpha set to: ", ALPHA_THRESHOLD, ") ---")
region_stats <- df_raw[, .(region_val = median(score, na.rm = TRUE)), by = .(region_id, type)]

# 1. Assumption Tests (using the Alpha Threshold)
shapiro_res <- region_stats %>% group_by(type) %>% shapiro_test(region_val)
normality <- all(shapiro_res$p > ALPHA_THRESHOLD)

levene_res <- region_stats %>% levene_test(region_val ~ type)
homogeneity <- levene_res$p > ALPHA_THRESHOLD

# 2. Test Selection and Post-hoc
if(normality && homogeneity) {
  test_used <- "ANOVA"
  post_hoc_res <- region_stats %>% tukey_hsd(region_val ~ type)
} else {
  test_used <- "Kruskal-Wallis"
  post_hoc_res <- region_stats %>% dunn_test(region_val ~ type, p.adjust.method = "bonferroni")
}

# 3. Customization of significance stars based on Alpha
# If p < ALPHA_THRESHOLD, we mark it as significant
post_hoc_res <- post_hoc_res %>%
  add_significance(cutpoints = c(0, 0.001, 0.01, ALPHA_THRESHOLD, 1),
                   symbols = c("***", "**", "*", "ns"))

write.csv(post_hoc_res, "LOH_posthoc_results.csv", row.names = FALSE)

# =============================
# Step 10: 📊 Plot B - Boxplot with geom_signif
# =============================
# Extracting comparisons and annotations for geom_signif
comp_list <- lapply(1:nrow(post_hoc_res), function(i) c(post_hoc_res$group1[i], post_hoc_res$group2[i]))
annot_list <- post_hoc_res$p.adj.signif

# Calculating bar heights to avoid overlap
max_val <- max(region_stats$region_val, na.rm = TRUE)
y_positions <- seq(max_val * 1.1, max_val * 1.5, length.out = length(comp_list))

pB <- ggplot(region_stats, aes(x = type, y = region_val, fill = type)) +
  geom_boxplot(width = 0.6, outlier.shape = 8, outlier.size = 1, alpha = 0.7) + 
  theme_classic() +
  scale_fill_manual(values = color_group) +
  geom_hline(yintercept = global_mean, linetype = "dashed", color = "red", size = 0.4, alpha = 0.5) +
  geom_signif(
    comparisons = comp_list,
    annotations = annot_list,
    y_position = y_positions,
    tip_length = 0.02,
    textsize = 4
  ) +
  
  scale_y_continuous(limits = c(0, max(y_positions) * 1.1)) +
  labs(x = NULL, y = "Median LOH Index per Region") +
  theme(legend.position = "none")
pB
# =============================
# Step 11: 🖼 Combine and Save
# =============================
final_plot <- plot_grid(pA, pB, ncol = 2, rel_heights = c(1.3, 1), labels = c("i", "ii"))
final_plot

ggsave("plotProfileTermination.pdf", pA, width = 4, height =3, units = "cm", dpi = 600)
ggsave("boxPlotTermination.pdf", pB, width = 10, height = 8, units = "cm", dpi = 600)

message("Process successfully completed!")
message("CSV File: LOH_posthoc_results.csv")
message("PDF File: LOH_Final_Profile_Stats.pdf")

