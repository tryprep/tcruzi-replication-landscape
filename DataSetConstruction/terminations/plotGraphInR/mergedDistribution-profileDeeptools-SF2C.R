# =============================
# Step 1: 📦 Install and load required packages
# =============================
packs <- c("ggplot2", "data.table", "R.utils", "patchwork")

missing <- packs[!(packs %in% installed.packages()[, "Package"])]
if (length(missing) > 0) {
  install.packages(missing, dependencies = TRUE)
}
suppressMessages(sapply(packs, require, character.only = TRUE))

# =============================
# Step 2: 📁 Set working directory
# =============================
setwd("~/Documents/DNAscent-Origins/dataAnalisys/Profile")

# =============================
# Step 3: 📂 Define input files
# =============================
matrix_files <- list(
  MoreFrequent = "matrix-moreFrequent-mergedTermination.gz",
  LessFrequent = "matrix-lessFrequent-mergedTermination.gz"
)

# =============================
# Step 4: 🔁 Define matrix processing function
# =============================
process_matrix_file <- function(file_path, bin_size = 10) {
  temp_file <- tempfile()
  gunzip(file_path, temp_file, remove = FALSE, overwrite = TRUE)
  
  df <- fread(temp_file, na.strings = c("nan", "NA"))
  unlink(temp_file)
  
  bins <- df[, 7:ncol(df), with = FALSE]
  mean_profile <- colMeans(bins, na.rm = TRUE)
  
  x_axis <- seq(from = 1, by = bin_size, length.out = length(mean_profile))
  
  return(data.frame(
    position = x_axis,
    mean_signal = mean_profile
  ))
}

# =============================
# Step 5: 📊 Define plotting function
# =============================
plot_profile <- function(plot_data, title_text, y_min = 0, y_max = 7, x_label_color = rgb(0, 0, 0, maxColorValue = 255)) {
  ggplot(plot_data, aes(x = position, y = mean_signal)) +
    geom_line(color = rgb(191,144,0, maxColorValue = 255), size = 1.5) +
    coord_cartesian(ylim = c(y_min, y_max)) +
    labs(
      x = NULL,
      y = "Count"
    ) +
    scale_x_continuous(
      breaks = c(min(plot_data$position), max(plot_data$position)),
      labels = c("Start", "End")
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.title = element_text(color ="black"),
      plot.title = element_text(hjust = 0.5)
    )
}

# =============================
# Step 6: 📈 Generate plots
# =============================
plot_data_list <- lapply(matrix_files, process_matrix_file)

plot_A <- plot_profile(
  plot_data_list$LessFrequent,
  "termination",
  x_label_color = rgb(26,76,166, maxColorValue = 255)
) +
  annotate("text", x = 0, y = Inf, label = "< 10", hjust = 0, vjust = 2,
           color = rgb(26,76,166, maxColorValue = 255), size = 5)

plot_B <- plot_profile(
  plot_data_list$MoreFrequent,
  "mini termination zones",
  x_label_color = rgb(200,35,0, maxColorValue = 255)
) +
  annotate("text", x = 0, y = Inf, label = ">= 10", hjust = 0, vjust = 2,
           color = rgb(200,35,0, maxColorValue = 255), size =5)

# =============================
# Step 7: 🧱 Combine and export plots
# =============================
final_plot <- plot_A + plot_B + plot_layout(ncol = 2)

# View plot
print(final_plot)

# Save plot to analysis folder
ggsave("terminationAndminiTerminantionzones.pdf", final_plot,  height = 9, width =16,dpi = 1200, units = "cm")
# Save plot to paper figure folder
setwd("~/Documents/paper_GenomeVariability/figures/Original/datasetConstrutionandClassification-D-NascentTermination")
ggsave("terminationAndminiTerminantionzones.pdf", final_plot,  height = 9, width =16,dpi = 1200, units = "cm")
