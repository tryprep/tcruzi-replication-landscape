# ==============================================================================
# R SCRIPT: GC/AT% Comparison (Custom Order + Pairwise Stats)
# ==============================================================================

# 0. Load Packages
library(tidyverse)
library(ggpubr)
library(rstatix)

# ==============================================================================
# STEP 1: Data Import
# ==============================================================================
file_path <- "~/Documents/DNAscent-Origins/dataAnalisys/CG-AT/all_origins_gc_at.csv"

if (file.exists(file_path)) {
  df <- read.csv(file_path)
  message("File loaded successfully.")
} else {
  stop("ERROR: File not found. Please check the path.")
}

# ==============================================================================
# STEP 2: Data Cleaning & Renaming (CRITICAL STEP)
# ==============================================================================

# 2.1 Define the raw names as they appear in your CSV file
raw_order <- c("singleOrigins", "consensusOrigins", 
               "singleTermination", "consensusTermination")

# 2.2 Define how you want them to appear in the plot (Labels)
pretty_labels <- c(
  "singleOrigins"        = "origins",
  "consensusOrigins"     = "IZD",
  "singleTermination"    = "termination",
  "consensusTermination" = "TZD"
)

# 2.3 Filter, Rename, and Set Factor Levels
df_clean <- df %>%
  # Keep only the groups we are interested in
  filter(origin_name %in% raw_order) %>%
  # CORREÇÃO AQUI: Forçando o uso do pacote dplyr explicitamente
  mutate(origin_name = dplyr::recode(origin_name, !!!pretty_labels)) %>%
  # Set the order (Levels) based on the new labels
  mutate(origin_name = factor(origin_name, levels = unname(pretty_labels)))

# Check if data is empty after filtering
if(nrow(df_clean) == 0) stop("Dataframe is empty. Check if CSV spelling matches 'raw_order'.")
# ==============================================================================
# STEP 3: Reshape Data (Wide to Long)
# ==============================================================================
df_long <- df_clean %>%
  pivot_longer(cols = "gc_percent", 
               names_to = "metric",
               values_to = "percent")

# ==============================================================================
# STEP 4: Normality Check
# ==============================================================================
normality_test <- df_long %>%
  group_by(metric, origin_name) %>%
  shapiro_test(percent)

print(normality_test)

# Decision logic: If p < 0.05 in ANY group, use Non-Parametric (Wilcoxon)
if (any(normality_test$p < 0.05)) {
  test_method <- "wilcox.test"
  test_name   <- "Wilcoxon (Non-Parametric)"
} else {
  test_method <- "t.test"
  test_name   <- "T-test (Parametric)"
}

# ==============================================================================
# STEP 5: Define Comparisons & Colors
# ==============================================================================

# 5.1 Comparisons List
# IMPORTANT: These names must match the "Pretty Labels" defined in Step 2.2
my_comparisons <- list(
  c("origins", "IZD"),
  c("origins", "termination"),
  c("origins", "TZD"),
  c("IZD", "termination"),
  c("IZD", "TZD"), 
  c("termination", "TZD")
)

# 5.2 Custom Colors
# Defining RGB colors for the "Pretty Labels"
origin_colors <- c( 
  "origins"          = rgb(116, 173, 209, maxColorValue = 255),  
  "termination"     = rgb(254, 224, 144, maxColorValue = 255),
  "IZD"              = rgb(244, 109, 67, maxColorValue = 255),
  "TZD"              = rgb(215, 48, 39, maxColorValue = 255)
)

# ==============================================================================
# STEP 6: Plotting
# ==============================================================================
plot_obj <- ggplot(df_long, aes(x = origin_name, y = percent, fill = origin_name)) +
  
  # Clean Boxplot (outliers hidden)
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.9) +
  
  # Mean Marker (White Diamond)
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  
  # Pairwise Statistics
  stat_compare_means(comparisons = my_comparisons,
                     method = test_method,
                     label = "p.signif", # ns, *, **, ***
                     tip.length = 0.01,
                     step.increase = 0.06, # Slightly increased to prevent overlap
                     hide.ns = FALSE, 
                     size= 3) +    # Change to TRUE if you want to hide 'ns'
  
  # Styling
  scale_fill_manual(values = origin_colors) + 
  labs(y = "Percentage (%)",
       x = NULL) + 
  
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.title.y = element_text(size = 12), 
        strip.text = element_text(size = 12, face = "bold"), 
        legend.position = "none") 

# ==============================================================================
# STEP 7: Export
# ==============================================================================
print(plot_obj)

# Define output directory
output_dir <- "~/Documents/paper_GenomeVariability/figures/Original/datasetConstrutionandClassification-D-NascentOrigins"


# Check if directory exists, if not, create it (optional safety)
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Save
ggsave(filename = file.path(output_dir, "GCcontent.pdf"), 
       plot = plot_obj, 
       width = 18, height = 14, units = "cm", dpi = 350)

message("Plot saved to: ", file.path(output_dir, "GCcontent.pdf"))
