# ==============================================================================
# SNP DENSITY ANALYSIS PIPELINE - LOCAL SERVER
# Description: This script processes genomic data to analyze SNP density across 
# different genomic compartments (Core, Disruptive, and GpDR).
# Input file: SNPs-compartment-Summary.csv
# ==============================================================================

# =======================================
# Step 1: 📦 Load required packages
# =======================================
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, dplyr, ggpubr)

# =======================================
# Step 2: 📁 Set working directory
# =======================================
setwd("~/Documents/DNAscent-Origins/dataAnalisys/Call-Variants/SNPs")

# =======================================
# Step 3: 📥 Import raw dataset 
# =======================================
data <- read.csv("SNPs-compartment-Summary.csv", header = TRUE, sep = "\t")
colnames(data) <- c("sample", "count", "length", "group", "SNPs_density")

# =======================================
# Step 4: 🔍 Data processing and normalization
# =======================================
# Function to filter and label compartments
process_comp <- function(data, total, without, single, consensus, label) {
  data %>%
    filter(group %in% c(total, without, single, consensus)) %>%
    mutate(group = case_when(
      group == total ~ "total_ref",
      group == without ~ "without\norigins",
      group == single ~ "origins",
      group == consensus ~ "IZD"
    )) %>%
    # Normalizing enrichment against the mean of the 'total_ref' group
    mutate(enrichment = SNPs_density / mean(SNPs_density[group == "total_ref"], na.rm = TRUE),
           compartment = label)
}

# Apply processing for each genomic compartment and combine data
df <- bind_rows(
  process_comp(data, "core", "core-Without-smAtlas", "origins-core", "miniInitiationZones-core", "Core"),
  process_comp(data, "disruptive", "disruptive-Without-smAtlas", "origins-disruptive", "miniInitiationZones-disruptive", "Disruptive"),
  process_comp(data, "GpDR", "GpDR-Without-smAtlas", "origins-GpDR", "miniInitiationZones-GpDR", "GpDR")
) %>%
  mutate(log2_enrichment = log2(enrichment)) %>%
  filter(is.finite(log2_enrichment), group != "total_ref")

# Set factor levels for consistent plotting order
df$group <- factor(df$group, levels = c("without\norigins", "origins", "IZD"))
df$compartment <- factor(df$compartment, levels = c("Core", "Disruptive", "GpDR"))
df$SNPsDensityPerkb <- as.numeric(df$SNPs_density * 1000)

# =======================================
# Step 5: 🎨 Define RGB color palette
# =======================================
group_colors <- c(
  "without\norigins" = rgb(0, 0, 255,   maxColorValue = 255), 
  "origins"          = rgb(0, 100, 255,  maxColorValue = 255),  
  "IZD"              = rgb(255, 0, 0,    maxColorValue = 255)
)

# =======================================
# Step 6: 🧪 Define pairwise comparisons
# =======================================
# Used for statistical mapping on the plot
my_comparisons <- list( 
  c("without\norigins", "origins"), 
  c("origins", "IZD"), 
  c("without\norigins", "IZD") 
)

# =======================================
# Step 7: 🖼️ Final Visualization
# =======================================
Plot_Final <- ggplot(df, aes(x = group, y = SNPsDensityPerkb, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
  facet_wrap(~compartment) +
  
  # Add Wilcoxon test results directly to the plot
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test", 
    label = "p.signif",
    step.increase = 0.12,
    tip.length = 0.01
  ) +
  
  theme_classic() +
  scale_fill_manual(values = group_colors) +
  labs(
    y = "SNP Density (SNPs/kb)",
    x = NULL
  ) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(color = "white", face = "bold", size = 12),
    # Formatting axis labels for better readability
    axis.text.x = element_text(size = 12, color = "black", lineheight = 0.9), 
    axis.title.y = element_text(size = 12), 
    axis.text.y = element_text(size = 12, color = "black", lineheight = 0.9)
  )

# Display the result
print(Plot_Final)

# =======================================
# Step 8: 💾 Save Output
# =======================================
# Save in the working directory
ggsave("SNPDensityInCompartment.pdf", Plot_Final, width = 16, height = 9, units = "cm", dpi = 300)

# Save in the paper figures directory
setwd("~/Documents/paper_GenomeVariability/figures/Original/SNPsDensity")
ggsave("SNPDensityInCompartment.pdf", Plot_Final, width = 18, height = 11, units = "cm", dpi = 300)