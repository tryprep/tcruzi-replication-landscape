# ==============================================================================
# Script Name: Genomic Compartment Origin Enrichment Visualization (Rose Plots)
# Description: This script calculates and visualizes the enrichment (fold change) 
#              of replication origins (single vs. consensus) across different 
#              T. cruzi genomic compartments (Core, Disruptive, GpDR) relative 
#              to a control dataset, using polar bar charts (rose plots).
#
# Required Input Files:
#   1. Control Dataset: "controlCompartment.txt" (Compartment vs. Peaks).
#   2. Experimental Datasets: 
#      - "NonSyncGenomeCompartment-lessFrequent.csv" (Single origins)
#      - "NonSyncGenomeCompartment-MoreFrequent.csv" (Consensus origins)
# ==============================================================================

# STEP 1: Load Dependencies
# Install any missing packages before loading
required_packages <- c("ggplot2", "dplyr", "rio", "scales", "cowplot")
missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]

if(length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}

# Load all required libraries silently
invisible(lapply(required_packages, require, character.only = TRUE))

# STEP 2: Environment Configuration
# Set your working directory
setwd("~/Documents/DNAscent-Origins/dataAnalisys/CompartmentAnnotation")

# STEP 3: Import and Standardize Data
# 3.1 Load Control
control <- rio::import("controlCompartment.txt")
colnames(control) <- c("Compartment", "Peaks")
control$Group <- "Control"

# 3.2 Load Experimental Datasets (Assigning group names directly)
nonsyncLess <- read.csv("NonSyncGenomeCompartment-lessFrequent.csv", header = FALSE, sep = "\t")
colnames(nonsyncLess) <- c("Compartment", "Peaks")
nonsyncLess$Group <- "single origins"

nonsyncMore <- read.csv("NonSyncGenomeCompartment-MoreFrequent.csv", header = FALSE, sep = "\t") 
colnames(nonsyncMore) <- c("Compartment", "Peaks")
nonsyncMore$Group <- "consensus origins"

# 3.3 Standardize Compartment names (Rename "Both" to "GpDR")
control$Compartment[control$Compartment == "Both"] <- "GpDR"
nonsyncLess$Compartment[nonsyncLess$Compartment == "Both"] <- "GpDR"
nonsyncMore$Compartment[nonsyncMore$Compartment == "Both"] <- "GpDR"

# STEP 4: Calculate Expected Frequencies (Control)
control_sums <- aggregate(Peaks ~ Compartment, data = control, FUN = sum)
control_sums$ExpectedFreq <- control_sums$Peaks / sum(control_sums$Peaks)

# STEP 5: Compute Observed Frequencies and Fold Change (Experimental)
experimental_all <- rbind(nonsyncLess, nonsyncMore)

# Sum Peaks per Group and Compartment
obs_sums <- aggregate(Peaks ~ Group + Compartment, data = experimental_all, FUN = sum)

# Sum total Peaks per Group for frequency denominator
group_totals <- aggregate(Peaks ~ Group, data = experimental_all, FUN = sum)
colnames(group_totals)[2] <- "TotalPeaks"

# Merge and calculate Fold Change
comparison_df <- merge(obs_sums, group_totals, by = "Group")
comparison_df$ObservedFreq <- comparison_df$Peaks / comparison_df$TotalPeaks
comparison_df <- merge(comparison_df, control_sums[, c("Compartment", "ExpectedFreq")], by = "Compartment")
comparison_df$FoldChange <- comparison_df$ObservedFreq / comparison_df$ExpectedFreq

# Factor levels for consistent plotting order
comparison_df$Compartment <- factor(comparison_df$Compartment, levels = c("Core", "Disruptive", "GpDR"))
comparison_df$Group <- factor(comparison_df$Group, levels = c("single origins", "consensus origins"))

# STEP 6: Define Custom Color Palettes
compartment_colors <- c(
  "Core"       = rgb(254, 196, 79, maxColorValue = 255),  # Light Gold (Mild)
  "Disruptive" = rgb(217, 95, 14, maxColorValue = 255),   # Intense Orange
  "GpDR"       = rgb(153, 52, 4, maxColorValue = 255)     # Burnt Orange/Brown (Very Intense)
)

group_text_colors <- c(
  "single origins"    = rgb(44, 123, 182, maxColorValue = 255),
  "consensus origins" = rgb(215, 25, 28, maxColorValue = 255)
)

# STEP 7: Define Plotting Function (Rose Plots)
create_rose_plot <- function(df_subset, group_name, text_color) {
  
  max_fc <- max(df_subset$FoldChange, na.rm = TRUE)
  y_limit <- max_fc * 1.35  # Buffer to prevent overlap with titles
  y_breaks <- seq(0, ceiling(max_fc), by = 1)
  
  ggplot(df_subset, aes(x = Compartment, y = FoldChange, fill = Compartment)) +
    geom_bar(stat = "identity", width = 1, color = "black") +
    coord_polar(start = 0) +
    scale_fill_manual(values = compartment_colors) +
    scale_y_continuous(limits = c(0, y_limit), breaks = y_breaks) +
    # Radial scale numbers
    geom_text(data = data.frame(y = y_breaks), 
              aes(x = 0.5, y = y, label = sprintf("%.1f", y)),
              inherit.aes = FALSE, size = 2, color = "black", vjust = -0.5) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid.major = element_line(color = "gray80", linewidth = 0.4), # Updated size to linewidth
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

# STEP 8: Generate Individual Plots
p1 <- create_rose_plot(subset(comparison_df, Group == "single origins"), 
                       "single origins", group_text_colors["single origins"])

p2 <- create_rose_plot(subset(comparison_df, Group == "consensus origins"), 
                       "consensus origins", group_text_colors["consensus origins"])

# STEP 9: Assemble and Save Final Figure
# 9.1 Extract the shared legend for compartments
dummy_for_legend <- ggplot(comparison_df, aes(x = Compartment, y = FoldChange, fill = Compartment)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Genome Compartments", values = compartment_colors) +
  theme_minimal() + 
  theme(legend.position = "bottom", 
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 10))

shared_legend <- cowplot::get_legend(dummy_for_legend)

# 9.2 Combine the two rose plots and assign facet labels I and II
combined_row <- cowplot::plot_grid(
  p1, p2, 
  ncol = 2, 
  labels = c("I", "II"), # Fixed missing labels argument
  label_size = 16, 
  label_fontface = "bold",
  label_x = 0.05, label_y = 0.95 # Position labels at top-left of each facet
)

# 9.3 Final assembly: Plots on top, legend on bottom
ROSES <- cowplot::plot_grid(combined_row, shared_legend, ncol = 1, rel_heights = c(1, 0.15))

# Display plot
print(ROSES)

# Save to PDF
ggsave("Figure_originsCompartmentAnnotated.pdf", 
       plot = ROSES, 
       height = 8, 
       width = 13, 
       units = "cm", 
       dpi = 1200)

message("Success: Compartment enrichment rose plots combined and saved.")
