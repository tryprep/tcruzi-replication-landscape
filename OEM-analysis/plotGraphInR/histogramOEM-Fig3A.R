# ==============================================================================
# Script Name: Genome-wide Origin Efficiency Metric (OEM) Distribution
# Description: This script visualizes the distribution of Origin Efficiency Metric 
#              (OEM) values across the genome. It generates a bar plot showing the 
#              binned frequencies of OEM values, overlaid with a LOESS smoothing 
#              curve to highlight the overall distribution trend.
#
# Required Input Files:
#   1. Binned Frequency Data: A CSV file named "data2-BinnedFrequency.csv" 
#      containing the distribution data (expected columns: 'Bin' and 'Frequency').
# ==============================================================================

# STEP 1: Load Dependencies
# Define required packages
packs <- c("ggplot2", "dplyr")

# Install missing packages
missing_packs <- packs[!packs %in% installed.packages()[, "Package"]]
if (length(missing_packs) > 0) {
  install.packages(missing_packs, dependencies = TRUE)
}

# Load all required packages silently
invisible(sapply(packs, require, character.only = TRUE))

# STEP 2: Environment Configuration
# Set your primary working directory
setwd("~/Documents/DNAscent-Origins/dataAnalisys/histogram")

# STEP 3: Load Data
# Read the binned frequency CSV file
df <- read.csv("data2-BinnedFrequency.csv")

# Uncomment the line below to inspect the data
# print(head(df))

# STEP 4: Data Visualization
# Create the distribution plot with a bar chart and a LOESS smoothing curve
plot_OEM <- ggplot(df, aes(x = Bin, y = Frequency)) +
  geom_bar(stat = "identity", fill = "darkgreen", width = 0.02) + 
  geom_smooth(aes(x = Bin, y = Frequency),
              method = "loess", se = FALSE,
              color = "blue", linewidth = 1.2, span = 0.8) +  # Updated 'size' to 'linewidth'
  labs(y = "Frequency (%)", x = "OEM Value", fill = "Experiment") +
  theme_classic() +
  coord_cartesian(ylim = c(0, 8)) +
  scale_y_continuous(breaks = seq(0, 8, by = 2),
                     expand = expansion(mult = c(0.005, 0.05))) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.25),
                     expand = expansion(mult = c(0.005, 0.05))) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 10, color = "black")
  )

# Display the plot in the R environment
print(plot_OEM)

# STEP 5: Save Plot
# 5.1 Save in the current analysis folder
ggsave("frequencyOEM.pdf", 
       plot = plot_OEM, 
       height = 8, 
       width = 9, 
       dpi = 1200, 
       units = "cm")
message("Saved: frequencyOEM.pdf in analysis folder.")

# 5.2 Save in the paper figure folder
setwd("~/Documents/paper_GenomeVariability/figures/Original/OEM")
ggsave("frequencyOEM.pdf", 
       plot = plot_OEM, 
       height = 8, 
       width = 9, 
       dpi = 1200, 
       units = "cm")
message("Saved: frequencyOEM.pdf in paper figures folder.")
