# ==============================================================================
# SNP DENSITY ANALYSIS PIPELINE - LOCAL SERVER
# Description: This script processes genomic data to analyze SNP density across 
# Orc1/Cdc6 binding sites.
# Input file: SNPs-InOriginsAndTerminationsAndOrcs-Summary.csv
# ==============================================================================

# =======================================
# Step 1: 📦 Load or install required packages
# =======================================
# Note: 'plyr' is excluded to avoid conflicts with 'dplyr'.
packs <- c("ggplot2", "rio", "dplyr", "ggsignif", 
           "rstatix", "RVAideMemoire", "car", "FSA", 
           "dunn.test", "patchwork")

# Check for missing packages and install them
missing_packs <- packs[!packs %in% installed.packages()[, "Package"]]
if (length(missing_packs) > 0) {
  install.packages(missing_packs, dependencies = TRUE)
}

# Load required packages
sapply(packs, require, character.only = TRUE)

# =======================================
# Step 2: 📁 Set working directory
# =======================================
# WARNING: Ensure this path exists on your local machine
setwd("~/Documents/DNAscent-Origins/dataAnalisys/Call-Variants/SNPs")

# =======================================
# Step 3: 📥 Import dataset
# =======================================
data <- read.csv("SNPs-InOriginsAndTerminationsAndOrcs-Summary.csv", 
                 header = TRUE, sep = "\t")

colnames(data) <- c("sample", "count", "length", "group", "SNPsDensity")

# =======================================
# Step 4: 🔍 Filter and Process Data
# =======================================
# Filter for Orc1/Cdc6 groups
Orcs_data <- data %>%
  filter(group %in% c("genome-2000-ALL-ORIGIN-FREE", "Orc1Cdc6"))

# Create the standardized column for Plotting (per kb)
# Calculated before rounding to maintain high precision
Orcs_data$SNPsPer1000bp <- Orcs_data$SNPsDensity * 1000

# Standardize group names (Recode for clearer plotting labels)
Orcs_data <- Orcs_data %>% 
  mutate(group = dplyr::recode(group, 
                               "genome-2000-ALL-ORIGIN-FREE" =  "Orc1/Cdc6 (-)", 
                               "Orc1Cdc6" = "Orc1/Cdc6 (+)"))

# Convert to data frame and define factor levels
df_Orcs <- as.data.frame(Orcs_data)
df_Orcs$group <- factor(df_Orcs$group,
                        levels = c("Orc1/Cdc6 (-)", "Orc1/Cdc6 (+)"))

# Inspect data structure
glimpse(df_Orcs)

# =======================================
# Step 5: 📊 Compute descriptive statistics
# =======================================
# Using dplyr for data summarization
mu_Orcs <- df_Orcs %>%
  group_by(group) %>%
  summarise(grp.mean = mean(SNPsPer1000bp, na.rm = TRUE),
            grp.median = median(SNPsPer1000bp, na.rm = TRUE))

print(mu_Orcs)

# =======================================
# Step 6: 📈 Shapiro-Wilk normality test
# =======================================
# Testing normality on the 'SNPsPer1000bp' column for each group
shapiro_results_Orcs <- sapply(levels(df_Orcs$group), function(g) {
  group_data <- df_Orcs$SNPsPer1000bp[df_Orcs$group == g]
  # Ensure sample size is sufficient for the test
  if (length(na.omit(group_data)) > 3) { 
    shapiro.test(group_data)$p.value
  } else {
    NA
  }
})

print(shapiro_results_Orcs)

# =======================================
# Step 7: 🧪 Choose and Run Statistical Test
# =======================================
# Select parametric or non-parametric test based on normality results (p > 0.05)
if (all(shapiro_results_Orcs[!is.na(shapiro_results_Orcs)] > 0.05)) {
  test_result_Orcs <- t.test(SNPsPer1000bp ~ group, data = df_Orcs)
  test_type_Orcs <- "t-test"
} else {
  test_result_Orcs <- wilcox.test(SNPsPer1000bp ~ group, data = df_Orcs)
  test_type_Orcs <- "Wilcoxon"
}

p_value_Orcs <- test_result_Orcs$p.value
cat("Test used:", test_type_Orcs, "| P-value:", p_value_Orcs, "\n")

# =======================================
# Step 8: ⭐ Convert p-value to significance symbol
# =======================================
p_to_symbol <- function(p) {
  if (is.na(p)) "ns"
  else if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else "ns"
}

signif_label_Orcs <- p_to_symbol(p_value_Orcs)

# =======================================
# Step 9: 🎨 Define colors for groups
# =======================================
# Colors are mapped directly to factor levels
colors_orcs <- c("Orc1/Cdc6 (-)" = rgb(44, 123, 182, maxColorValue = 255),
                 "Orc1/Cdc6 (+)"   = rgb(215, 25, 28, maxColorValue = 255))

# =======================================
# Step 10: 🖼️ Create Box plot
# =======================================
p <- ggplot(df_Orcs, aes(x = group, y = SNPsPer1000bp, fill = group)) +
  geom_boxplot(width = 0.8,         
               lwd = 0.5,           # Use 'lwd' for cross-version compatibility
               outlier.shape = NA, 
               alpha = 1) +
  theme_classic() +
  labs(y = "SNP density (SNPs/kb)", x = NULL) + 
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, angle = 0, vjust = 0.7),
        axis.title.y = element_text(vjust = 4, size = 12),
        legend.position = "none",
        plot.margin = margin(t = 10, r = 10, b = 15, l = 20),
        panel.spacing = unit(1, "lines")) +
  scale_fill_manual(values = colors_orcs) +
  # Add statistical significance annotation
  geom_signif(comparisons = list(c("Orc1/Cdc6 (-)", "Orc1/Cdc6 (+)")),
              annotations = signif_label_Orcs,
              y_position = 14,      # Adjust based on data range
              tip_length = 0.02,
              size = 0.5,           # Line thickness       
              textsize = 5.0,
              vjust = -0.2) +       # Adjust text vertical alignment
  ylim(11, 15) 

# =======================================
# Step 11: 📊 Display plot
# =======================================
print(p)

# =======================================
# Step 12: 💾 Export plot to PDF
# =======================================
output_path_1 <- "~/Documents/DNAscent-Origins/dataAnalisys/Call-Variants/SNPs"
output_path_2 <- "~/Documents/paper_GenomeVariability/figures/Original/SNPsDensity"

# Save to the first directory if it exists
if(dir.exists(output_path_1)) {
  ggsave(filename = file.path(output_path_1, "SNPsInOrc1Cdc6BindingSites.pdf"), 
         plot = p, height = 10, width = 12, dpi = 1200, units = "cm")
}

# Save to the second directory (Manuscript Figures) if it exists
if(dir.exists(output_path_2)) {
  ggsave(filename = file.path(output_path_2, "SNPsInOrc1Cdc6BindingSites.pdf"), 
         plot = p, height = 9, width = 10, dpi = 1200, units = "cm")
}
