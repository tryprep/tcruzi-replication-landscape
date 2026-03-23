# ==============================================================================
# ANALYSIS OVERVIEW:
# This plot represents the SNP density profile along a region surrounding the 
# origin of replication or termination, oriented according to the transcription 
# direction of the genes within the analyzed coordinates.
# 
# This analysis was performed on the local server, following the variant calling 
# and SNP analysis pipeline previously mentioned.
# ==============================================================================

# ==============================================================================
# Step 1: 📦 Load or install required packages
# ==============================================================================
packs <- c("ggplot2", "dplyr", "ggsignif", "tidybayes", "tidyr")
# Note: plyr and FSA were removed if not strictly used to avoid namespace conflicts
lapply(packs, require, character.only = TRUE)

# ==============================================================================
# Step 2: Set working directory
# ==============================================================================
setwd("~/Documents/DNAscent-Origins/dataAnalisys/Call-Variants/SNPs")

# ==============================================================================
# Step 3: 📂 Read input data and select group
# ==============================================================================
# NOTE: The "ORC_regions_table_2000bp.tsv" file contains SNP density information 
# for all analyzed groups (e.g., origins, IZD, terminations, TZD).
# To plot a specific group, simply change the 'selected_group' variable below.

selected_group <- "miniInitiationZones" # e.g., "origins", "IZD", "terminations", "TZD"

rawData <- read.csv("ORC_regions_table_2000bp.tsv", header = TRUE, sep = "\t", quote = "\t")

dataset_filtered <- rawData %>% 
  filter(
    Group == selected_group , 
    Strand %in% c("-", "+") , 
    SNPs_density > 0 )

# ==============================================================================
# Step 4: 🧹 Clean and process data 
# ==============================================================================
data <- as.data.frame(dataset_filtered)
# Normalize density to per kb and clean Strand strings
data$SNPs_density <- data$SNPs_density * 1000
data$Strand <- gsub("[[:space:]]", "", as.character(data$Strand))

data <- data %>%
  mutate(
    # Extract numeric value from the region string (e.g., "upstream_2000" -> 2000)
    dist_val = suppressWarnings(as.numeric(gsub("[^0-9]", "", region))),
    # Define actual position: Upstream is negative, Downstream is positive, Site is 0
    region_pos = case_when(
      grepl("upstream", region) ~ -dist_val,
      grepl("downstream", region) ~ dist_val,
      region == "site" ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  # Keep only the window of 12kb around the origin (from -6000 to 6000)
  filter(!is.na(region_pos), abs(region_pos) <= 6000) %>%
  # Maintain only valid strands and regions with minimum length
  filter(Strand %in% c("+", "-"), length >= 2000)

# Save the processed dataframe
write.csv(data, "processed_snp_density_data.csv", row.names = FALSE)

# ==============================================================================
# Step 5: 🧪 Create Stratified Bootstrap dataframe
# ==============================================================================
n_reps <- 15000
boot_results <- data %>%
  filter(!is.na(SNPs_density)) %>% 
  group_by(region_pos, Strand) %>%
  reframe({
    
    s_list <- split(SNPs_density, sample)
    
    b_dist <- replicate(n_reps, {
      sample_means <- sapply(s_list, function(x) {
        if(length(x) > 0) mean(sample(x, replace = TRUE), na.rm = TRUE) else NA
      })
      mean(sample_means, na.rm = TRUE)
    })
    
    if(all(is.na(b_dist)) || var(b_dist, na.rm = TRUE) == 0) {
      m <- median(b_dist, na.rm = TRUE)
      h_low <- m
      h_high <- m
    } else {
      interval <- tidybayes::hdi(b_dist, .width = 0.94 , na.rm = TRUE)
      h_low  <- as.numeric(min(interval[,1]))
      h_high <- as.numeric(max(interval[,2]))
    }
    
    tibble(
      SNPs_density = median(b_dist, na.rm = TRUE),
      hdi_low = h_low,
      hdi_high = h_high
    )
  }) %>%
  filter(!is.na(SNPs_density))

# ==============================================================================
# Step 6: 📈 Plot Function 
# ==============================================================================
############################################################################################################################
plot_snp_wave_smooth <- function(df, y_limits = c(11.5, 13.5), smoothness = 0.3) {
  
  # 1. PRE-PROCESSING: Manual smoothing of the 3 columns
  # This ensures that the ribbon (HDI) is drawn perfectly
  df_smooth <- df %>%
    group_by(Strand) %>%
    arrange(region_pos) %>%
    do({
      dat <- .
      # Create a denser grid of points to make the curve "smooth"
      x_grid <- seq(min(dat$region_pos), max(dat$region_pos), length.out = 300)
      
      # Apply loess individually for each metric
      # The span parameter controls the "smoothness"
      fit_y    <- loess(SNPs_density ~ region_pos, data = dat, span = smoothness)
      fit_low  <- loess(hdi_low ~ region_pos, data = dat, span = smoothness)
      fit_high <- loess(hdi_high ~ region_pos, data = dat, span = smoothness)
      
      data.frame(
        region_pos = x_grid,
        SNPs_density = predict(fit_y, x_grid),
        hdi_low = predict(fit_low, x_grid),
        hdi_high = predict(fit_high, x_grid)
      )
    }) %>%
    ungroup()
  
  # 2. PLOTTING
  strand_colors <- c("+" = "blue", "-" = "red")
  strand_fills  <- c("+" = "#6495ED", "-" = "#FFA07A")
  y_annot <- y_limits[2]
  
  plt <- ggplot(df_smooth, aes(x = region_pos, group = Strand, color = Strand, fill = Strand)) +
    
    # Ribbon (HDI) - This will now appear since the data is ready
    geom_ribbon(aes(ymin = hdi_low, ymax = hdi_high), alpha = 0.3, color = NA) +
    
    # Mean line
    geom_line(aes(y = SNPs_density), linewidth = 1.1) +
    
    # References (Site)
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
    # Upstream Boundary
    geom_vline(xintercept = -6000, linetype = "dashed", color = "grey70") +
    annotate("text", x = -4000, y = 11.5, label = "Upstream", angle = 0, vjust = -0.5, color = "black", size = 3) +
    
    # Downstream Boundary
    geom_vline(xintercept = 6000, linetype = "dashed", color = "grey70") +
    annotate("text", x = 4000, y = 11.5, label = "Downstream", angle = 0, vjust = -0.5, color = "black", size = 3) +
    
    # Axes and Scales
    scale_x_continuous(
      breaks = seq(-6000, 6000, by = 2000),
      labels = function(x) ifelse(x == 0, "Site", paste0(x / 1000))
    ) +
    scale_color_manual(values = strand_colors) +
    scale_fill_manual(values = strand_fills) +
    
    coord_cartesian(ylim = y_limits, expand = FALSE) +
    labs(
      x = "Distance from origin site (kb)", 
      y = "SNP Density per kb",
      color = "Strand", fill = "Strand"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey95"), 
      axis.text = element_text(size = 12)
    )
  
  return(plt)
}

# ==============================================================================
# Step 7: 🚀 Run the function
# ==============================================================================
# Use y_limits = c(10, 13) based on your data summary
final_plot <- plot_snp_wave_smooth(boot_results, y_limits = c(10, 40), smoothness = 0.3)
print(final_plot)

# ==============================================================================
# Step 8: 🚀 Save to plot profile (Dynamic Filename)
# ==============================================================================
# The filename will automatically update based on the chosen 'selected_group'
output_filename <- paste0("SNPsDensityPerStrands_", selected_group, ".pdf")
ggsave(output_filename, plot = final_plot, width = 9, height = 7, units = "cm", dpi = 1200)

# ==============================================================================
# Step 9: 🚀 Plot zoom 
# ==============================================================================
# Creating the zoom dataframe with specific points
zoom <- boot_results %>%
  filter(region_pos %in% c(-2000, 0, 2000))

# View the result to check
print(zoom)

plot_snp_wave_zoom <- function(df, y_limits = c(11.5, 13.5), smoothness = 0.8) {
  
  # 1. Manual smoothing (necessary for loess to work with few points)
  df_smooth <- df %>%
    group_by(Strand) %>%
    arrange(region_pos) %>%
    do({
      dat <- .
      x_grid <- seq(min(dat$region_pos), max(dat$region_pos), length.out = 200)
      
      fit_y    <- loess(SNPs_density ~ region_pos, data = dat, span = smoothness)
      fit_low  <- loess(hdi_low ~ region_pos, data = dat, span = smoothness)
      fit_high <- loess(hdi_high ~ region_pos, data = dat, span = smoothness)
      
      data.frame(
        region_pos = x_grid,
        SNPs_density = predict(fit_y, x_grid),
        hdi_low = predict(fit_low, x_grid),
        hdi_high = predict(fit_high, x_grid)
      )
    }) %>%
    ungroup()
  
  strand_colors <- c("+" = "blue", "-" = "red")
  strand_fills  <- c("+" = "#6495ED", "-" = "#FFA07A")
  
  plt <- ggplot(df_smooth, aes(x = region_pos, group = Strand, color = Strand, fill = Strand)) +
    
    # Ribbon and Smooth Line
    geom_ribbon(aes(ymin = hdi_low, ymax = hdi_high), alpha = 0.3, color = NA) +
    geom_line(aes(y = SNPs_density), linewidth = 0.7) +
    
    # --- Requested Annotations ---
    # Central Reference (Site)
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.5, linewidth = 0.3) +
    
    # Upstream
    geom_vline(xintercept = -2000, linetype = "dashed", color = "grey70", linewidth = 0.3) +
    
    # Downstream
    geom_vline(xintercept = 2000, linetype = "dashed", color = "grey70", linewidth = 0.3) +
    
    # --- Axis and Layout Configurations ---
    scale_x_continuous(
      breaks = c(-2000, 0, 2000),
      labels = c("-2", "Site", "2")
    ) +
    scale_color_manual(values = strand_colors) +
    scale_fill_manual(values = strand_fills) +
    
    coord_cartesian(ylim = y_limits, expand = FALSE) +
    labs( x = "Distance from \n origin site (kb)", 
          y = "SNP Density \n per kb") +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey95"),
      axis.title = element_text(size = 5),
      axis.text = element_text(size = 4),
      plot.margin = margin(5, 10, 5, 10) # Margins so the text is not cut off
    )
  
  return(plt)
}

# To view now:
final_plot_zoom <- plot_snp_wave_zoom(zoom, y_limits = c(10,30))
print(final_plot_zoom)

# Save to PDF (Dynamic Filename)
zoom_output_filename <- paste0("zoom-single_", selected_group, ".pdf")
ggsave(zoom_output_filename, 
       plot = final_plot_zoom, 
       width = 4, 
       height = 3, 
       units = "cm",
       device = "pdf", 
       bg = "transparent")
