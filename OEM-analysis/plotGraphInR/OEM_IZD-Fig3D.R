# ==============================================================================
# Script Name: Initiation Zone Domain (IZD) OEM Profile Analysis
# Description: This script computes and visualizes the Origin Efficiency Metric 
#              (OEM) profile across merged regions that define Initiation Zone 
#              Domains (IZDs). It handles edge cases (padding out-of-bounds 
#              regions with NAs), dynamically bins the data, and generates both 
#              a "spaghetti plot" of individual profiles and a zoomed-in mean 
#              metaplot.
#
# Required Input Files:
#   1. Merged BED File: "consensus-merged.bed" containing the IZD coordinates.
#   2. OEM bedGraph: "OEM.bedgraph" containing genome-wide OEM scores.
# ==============================================================================

# STEP 1: Load Dependencies
required_packages <- c("data.table", "ggplot2", "cowplot")
bioc_packages <- c("GenomicRanges", "rtracklayer", "IRanges")

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

# STEP 2: Environment Configuration & Parameters
SUBDIR <- "OEM_Analisys/"
base_dir <- file.path("~/Documents/DNAscent-Origins/dataAnalisys/Profile", SUBDIR)

# Ensure the directory exists before attempting to set the working directory
if (!dir.exists(base_dir)) {
  dir.create(base_dir, recursive = TRUE)
}
setwd(base_dir)

# --- FILES ---
my_bed_file   <- "consensus-merged.bed"        
bedgraph_file <- "OEM.bedgraph"        

# --- PARAMETERS ---
bin_size_bp  <- 20    # Bin size of 20bp
plot_radius  <- 20000 # 20kb to each side (Total window 40kb)

# STEP 3: Load & Process Data
message(">> Loading data...")
bg <- fread(bedgraph_file, header = FALSE)
setnames(bg, c("chr", "start", "end", "score"))

# Note: BedGraph is 0-based start, 1-based end. GRanges is 1-based.
bg.gr <- GRanges(seqnames = bg$chr, 
                 ranges = IRanges(start = bg$start + 1, end = bg$end), 
                 score = bg$score)

bed <- fread(my_bed_file, header = FALSE)
setnames(bed, 1:3, c("chr", "start", "end"))
bed_raw_gr <- GRanges(seqnames = bed$chr, ranges = IRanges(start = bed$start + 1, end = bed$end))

# Prepare Coverage
cvg <- coverage(bg.gr, weight = "score")

# Sync Chromosomes
common_chr <- intersect(seqlevels(bed_raw_gr), names(cvg))
if (length(common_chr) == 0) stop("ERROR: No common chromosomes found.")
bed_raw_gr <- keepSeqlevels(bed_raw_gr, common_chr, pruning.mode = "coarse")

# --- CENTER REGIONS ---
message(">> Centering windows...")
centers <- resize(bed_raw_gr, width = 1, fix = "center")

# Define total target width
target_width <- plot_radius * 2
windows <- resize(centers, width = target_width, fix = "center")

# Inform windows about chromosome lengths (crucial for boundary checks)
seqlengths(windows) <- elementNROWS(cvg)[seqlevels(windows)]

# NOTE: We do NOT use trim() here. We keep the "ideal" coordinates 
# even if they extend beyond the chromosome boundaries.
message(paste(">> Total regions for plotting (including edges):", length(windows)))

# STEP 4: Extraction, Padding & Binning
message(">> Extracting, Padding (NA), and Binning...")

# Function to bin vector (handles NA by ignoring them in mean)
bin_vector_mean <- function(vec, bsize) {
  if(length(vec) == 0) return(numeric(0))
  grouping <- rep(1:ceiling(length(vec)/bsize), each = bsize)[1:length(vec)]
  tapply(vec, grouping, mean, na.rm = TRUE)
}

all_binned_profiles <- list()

for (chrom in common_chr) {
  chr_wins <- windows[seqnames(windows) == chrom]
  if (length(chr_wins) == 0) next
  
  chr_cvg <- cvg[[chrom]]
  chr_len <- length(chr_cvg)
  
  # 1. Identify "Ideal" coordinates (requested window)
  ideal_starts <- start(chr_wins)
  ideal_ends   <- end(chr_wins)
  
  # 2. Identify "Valid" coordinates (clamped to chromosome limits)
  valid_starts <- pmax(1, ideal_starts)
  valid_ends   <- pmin(chr_len, ideal_ends)
  
  # 3. Extract data only from valid ranges
  # 'Views' will fail if coordinates are out of bounds, so we use the clamped ranges
  valid_ranges <- IRanges(start = valid_starts, end = valid_ends)
  chr_views <- Views(chr_cvg, valid_ranges)
  raw_data_list <- lapply(chr_views, as.numeric)
  
  # 4. Apply Padding (Fill missing edges with NA)
  # mapply processes each region individually
  padded_data_list <- mapply(function(d, i_start, i_end, v_start, v_end) {
    
    pad_left  <- v_start - i_start # How many bases missing on the left?
    pad_right <- i_end - v_end     # How many bases missing on the right?
    
    # Reconstruct: [NAs Left] + [Real Data] + [NAs Right]
    c(rep(NA, pad_left), d, rep(NA, pad_right))
    
  }, raw_data_list, ideal_starts, ideal_ends, valid_starts, valid_ends, SIMPLIFY = FALSE)
  
  # 5. Binning the padded vectors
  chr_binned <- lapply(padded_data_list, bin_vector_mean, bsize = bin_size_bp)
  
  all_binned_profiles <- c(all_binned_profiles, chr_binned)
}

# Consistency Check: Ensure all rows have the exact same number of bins
expected_bins <- ceiling(target_width / bin_size_bp)
lens <- sapply(all_binned_profiles, length)

if(any(lens != expected_bins)) {
  warning("Some rows have different lengths due to rounding. Trimming to match expected size.")
  all_binned_profiles <- lapply(all_binned_profiles, function(x) x[1:expected_bins])
}

mat <- do.call(rbind, all_binned_profiles)

if (is.null(mat) || nrow(mat) == 0) stop("ERROR: Empty matrix.")

# STEP 5: Data Reshaping
message(">> Formatting data structure...")

df_heat <- as.data.table(mat)

# --- Fix for 'Pattern not found' in melt ---
# Explicitly renaming columns to V1, V2, V3... to ensure melt works reliably
message(">> Renaming columns for compatibility...")
new_col_names <- paste0("V", 1:ncol(df_heat))
setnames(df_heat, new_col_names)

# Add ID
df_heat[, region_id := 1:.N]

# Melt (Safe mode)
df_long <- melt(df_heat, 
                id.vars = "region_id", 
                measure.vars = patterns("^V"), # Matches V1, V2...
                variable.name = "bin_label", 
                value.name = "signal")

# Convert label V1... to number, then to Distance
df_long[, bin_idx := as.numeric(gsub("V", "", bin_label))]

# Calculate X Axis (Distance from Center)
# Formula: (Bin Index * Bin Size) - Radius - (Half bin for centering)
df_long[, dist_bp := (bin_idx * bin_size_bp) - plot_radius - (bin_size_bp / 2)]

# Global Mean Calculation
# na.rm = TRUE is critical here because padded regions contain NAs
df_mean <- df_long[, .(mean_signal = mean(signal, na.rm = TRUE)), by = dist_bp]

# STEP 6: Visualization
message(">> Generating plots...")

n_lines <- nrow(mat)
# Dynamic transparency adjustment based on number of regions
my_alpha <- if(n_lines > 5000) 0.02 else if(n_lines > 1000) 0.05 else 0.3

# Plot 1: Centered Spaghetti Plot
p <- ggplot() +
  # Individual lines (Gray)
  # Gaps (NA) will break the lines, which is the correct behavior for edges
  geom_line(data = df_long, 
            aes(x = dist_bp, y = signal, group = region_id), 
            color = "gray20", 
            alpha = my_alpha, 
            linewidth = 0.1) +
  
  # Mean line (Red)
  geom_line(data = df_mean,
            aes(x = dist_bp, y = mean_signal),
            color = "firebrick",
            linewidth = 1.2) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.6) +
  geom_vline(xintercept = 2000, linetype = "dashed", color = "blue", alpha = 0.6) +
  geom_vline(xintercept = -2000, linetype = "dashed", color = "blue", alpha = 0.6) +
  
  scale_x_continuous(name = "Distance from Center (bp)",
                     breaks = seq(-plot_radius, plot_radius, by = 5000)) + 
  
  labs(y = "OEM Signal") +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.title = element_text(face = "bold"), 
    axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12), 
    axis.title = element_text(size = 12)
  ) 

print(p)

# Plot 2: Zoomed Mean Plot
z <- ggplot() +
  # Keep only the mean line (Red)
  geom_line(data = df_mean,
            aes(x = dist_bp, y = mean_signal),
            color = "firebrick",
            linewidth = 1.2) +
  
  # Vertical center lines
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.6) +
  geom_vline(xintercept = 2000, linetype = "dashed", color = "blue", alpha = 0.6) +
  geom_vline(xintercept = -2000, linetype = "dashed", color = "blue", alpha = 0.6) +
  
  # Zoom on the Y-axis (without dropping data)
  coord_cartesian(ylim = c(-0.25, 0.25)) + 
  
  # X-axis configuration
  scale_x_continuous(name = "Distance from Center (bp)",
                     breaks = seq(-plot_radius, plot_radius, by = 5000)) +
  
  labs(y = "OEM Signal") +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
    
    # Reducing all legends and texts to size 6
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.key.size = unit(0.2, "cm"),
    
    # Axis texts reduced for the zoom plot
    axis.text.x = element_text(angle = 60, hjust = 1, size = 6, color = "black"),
    axis.text.y = element_text(size = 6, color = "black"),
    axis.title = element_text(size = 7, face = "bold"),
    
    # Removing background grid for a cleaner look
    panel.grid.minor = element_blank()
  )

print(z)

# STEP 7: Save Plots
out_file <- paste0("Centered_Spaghetti_Padded_", tools::file_path_sans_ext(my_bed_file), ".pdf")
ggsave(out_file, plot = p, width = 20, height = 9, units = "cm", dpi = 1200)

ggsave("zoom-OEM_Centered_Spaghetti_Padded.pdf", plot = z, width = 5, height = 3, units = "cm", dpi = 1200)

message(">> Success: Plots saved as ", out_file, " and zoom-OEM_Centered_Spaghetti_Padded.pdf")