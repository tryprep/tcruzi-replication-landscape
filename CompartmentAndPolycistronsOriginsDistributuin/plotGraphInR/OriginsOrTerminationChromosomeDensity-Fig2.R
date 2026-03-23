# ==============================================================================
# Script Name: Chromosome Origin Density Visualization
# Description: This script calculates and visualizes the density (bp/bp) of 
#              replication origins or terminations across chromosomes. It computes 
#              the total length of origins/terminations per chromosome and divides 
#              it by the respective chromosome size to generate density bar plots.
#
# Required Input Files:
#   1. BED Files: Containing the genomic coordinates of the replication origins 
#      or terminations (e.g., "control.bed", "singleOrigins.bed"). Must have at 
#      least 3 columns: Chromosome, Start, and End.
#   2. Chromosome Sizes File: A text file named "genome-chromoSizes.txt" 
#      containing two tab-separated columns: Chromosome Name and Size.
# ==============================================================================

# STEP 1: Load Dependencies
library(tidyverse)
library(stringr)

# STEP 2: Environment Configuration
# Set your working directory
setwd("~/Documents/DNAscent-Origins/dataAnalisys/ChromossomeDistribution")

# Load chromosome sizes (Expected format: Name \t Size)
chrom_sizes <- read.table("genome-chromoSizes.txt", 
                          header = FALSE, 
                          col.names = c("Chromosome", "ChrSize"),
                          stringsAsFactors = FALSE)

# STEP 3: Define Group Mapping and Labels
# Using \n for line breaks in facet labels to keep the plot compact
group_mapping <- c(
  "control.bed"               = "Without \n Origins",
  "singleOrigins.bed"         = "Single \n Origins",
  "consensusOrigins.bed"      = "Consensus \n Origins",
  "singleTerminations.bed"    = "Single \n Terminations",
  "consensusTerminations.bed" = "Consensus \n Terminations" 
)

files <- names(group_mapping)

# STEP 4: Define Data Processing Function
process_origins <- function(file_path) {
  # Check if file exists to prevent crash
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  # Read BED file (assuming at least 3 columns: Chrom, Start, End)
  data <- read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  data %>%
    dplyr::select(V1, V2, V3) %>%
    dplyr::rename(Chromosome = V1, Start = V2, End = V3) %>%
    dplyr::mutate(
      Length = End - Start,
      GroupFile = basename(file_path)
    )
}

# STEP 5: Import and Merge Data
# map_df automatically binds the rows of the 5 files
all_origins <- map_df(files, process_origins)

# STEP 6: Sorting Logic and Density Calculation
# Extract numeric ID from "TcChr41-S" to sort 1 to 41 correctly
chrom_order <- chrom_sizes %>%
  dplyr::mutate(ChrNum = as.numeric(str_extract(Chromosome, "\\d+"))) %>%
  dplyr::arrange(ChrNum) %>%
  dplyr::pull(Chromosome)

density_data <- all_origins %>%
  dplyr::group_by(Chromosome, GroupFile) %>%
  dplyr::summarise(Total_Origin_Length = sum(Length), .groups = "drop") %>%
  dplyr::inner_join(chrom_sizes, by = "Chromosome") %>%
  dplyr::mutate(
    Density = Total_Origin_Length / ChrSize,
    # Assign labels and lock order for plotting
    GroupName = factor(group_mapping[GroupFile], levels = unname(group_mapping)),
    # Lock numeric order for Chromosomes
    Chromosome = factor(Chromosome, levels = chrom_order)
  )

# STEP 7: Define Visualization Parameters
# Mapping specific hex/rgb colors to the group labels
origin_colors <- c( 
  "Without \n Origins"          = rgb(0, 0, 255,   maxColorValue = 255), 
  "Single \n Origins"           = rgb(0, 100, 255, maxColorValue = 255),  
  "Consensus \n Origins"        = rgb(255, 0, 0,   maxColorValue = 255), 
  "Single \n Terminations"      = rgb(100, 100, 255, maxColorValue = 255),
  "Consensus \n Terminations"   = rgb(255, 25, 0,  maxColorValue = 255)
)

# STEP 8 & 9: Create and Save Individual Plots
# Get the list of unique groups respecting the factor order
unique_groups <- levels(density_data$GroupName)

for (group in unique_groups) {
  plot_data <- density_data %>% 
    dplyr::filter(GroupName == group)
  
  plot_density <- ggplot(plot_data, aes(x = Chromosome, y = Density, fill = GroupName)) +
    geom_bar(stat = "identity", width = 0.7) +
    theme_minimal() +
    scale_fill_manual(values = origin_colors) + 
    labs(
      x = NULL,
      y = "Density (bp / bp)"
    ) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      panel.grid.minor = element_blank(),
      legend.position = "none" # Since there is only one group per plot, the legend is unnecessary
    ) 
  # + coord_flip()
  
  print(plot_density)
  
  # Replace line breaks and spaces for a clean filename
  safe_filename <- gsub(" \n | ", "_", group)
  file_name <- paste0("Origin_Density_", safe_filename, ".pdf")
  
  ggsave(filename = file_name, 
         plot = plot_density, 
         width = 10, 
         height = 8, # Reduced height since there are no stacked facets
         dpi = 1200, 
         units = "cm")
  
  message(paste("Saved:", file_name))
}

message("Success: All individual plots generated and saved.")