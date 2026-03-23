# ==============================================================================
# SNP DENSITY ANALYSIS PIPELINE - LOCAL SERVER
# Description: This script processes genomic data to analyze SNP density across 
# different Orc1/Cdc6 binding sites (Active, Inactive, and Origin-free).
# Input file: SNPs-InOriginsAndTerminationsAndOrcs-Summary.csv
# ==============================================================================

# =======================================
# Step 1: 📦 Load or install required packages
# =======================================
packs <- c("ggplot2", "rio", "dplyr", "ggsignif", 
           "rstatix", "RVAideMemoire", "car", "FSA", 
           "dunn.test", "patchwork")

missing_packs <- packs[!packs %in% installed.packages()[, "Package"]]
if (length(missing_packs) > 0) {
  install.packages(missing_packs, dependencies = TRUE)
}

sapply(packs, require, character.only = TRUE)

# =======================================
# Step 2: 📁 Set working directory
# =======================================
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
Orcs_data <- data %>%
  filter(group %in% c("genome-2000-ALL-ORIGIN-FREE", 
                      "Orc1Cdc6-onActiveOrigins", 
                      "Orc1Cdc6-onInactiveOrigins"))

# Create the standardized column (per kb)
Orcs_data$SNPsPer1000bp <- Orcs_data$SNPsDensity * 1000

# Recode for better plotting labels
Orcs_data <- Orcs_data %>% 
  mutate(group = dplyr::recode(group, 
                               "genome-2000-ALL-ORIGIN-FREE" = "Orc1/Cdc6 (-)",
                               "Orc1Cdc6-onInactiveOrigins"  = "Orc1/Cdc6 (+) \n (-) ",
                               "Orc1Cdc6-onActiveOrigins"    = "Orc1/Cdc6 (+) \n (+)"))

# Convert to factor and define order
df_Orcs <- as.data.frame(Orcs_data)
df_Orcs$group <- factor(df_Orcs$group,
                        levels = c("Orc1/Cdc6 (-)", "Orc1/Cdc6 (+) \n (-) ","Orc1/Cdc6 (+) \n (+)"))

glimpse(df_Orcs)

# =======================================
# Step 5: 📊 Compute descriptive statistics
# =======================================
mu_Orcs <- df_Orcs %>%
  group_by(group) %>%
  summarise(grp.mean = mean(SNPsPer1000bp, na.rm = TRUE),
            grp.median = median(SNPsPer1000bp, na.rm = TRUE),
            n = n())

print(mu_Orcs)

# =======================================
# Step 6: 📈 Shapiro-Wilk normality test
# =======================================
shapiro_results_Orcs <- df_Orcs %>%
  group_by(group) %>%
  shapiro_test(SNPsPer1000bp)

print(shapiro_results_Orcs)

# =======================================
# Step 7: 🧪 Statistical Test (3 groups)
# =======================================
# Check normality (p > 0.05 for ALL groups)
is_normal <- all(shapiro_results_Orcs$p > 0.05)

if (is_normal) {
  res_main <- aov(SNPsPer1000bp ~ group, data = df_Orcs)
  # Native R alternative:
  post_hoc_raw <- TukeyHSD(res_main)
  print(post_hoc_raw)
  
  # Ensure compatibility for geom_signif plotting:
  post_hoc <- as.data.frame(post_hoc_raw$group) 
  test_type <- "ANOVA" 
} else {
  # Non-parametric: Kruskal-Wallis + Dunn's test
  res_main <- kruskal.test(SNPsPer1000bp ~ group, data = df_Orcs)
  post_hoc <- df_Orcs %>% dunn_test(SNPsPer1000bp ~ group, p.adjust.method = "bonferroni")
  test_type <- "Kruskal-Wallis"
}

print(test_type)
print(post_hoc)

# =======================================
# Step 8: ⭐ Define Comparisons for Plotting
# =======================================
# Define which pairs you want to compare in the plot
my_comparisons <- list(
  c("Orc1/Cdc6 (-)", "Orc1/Cdc6 (+) \n (-) "),
  c("Orc1/Cdc6 (+) \n (-) ", "Orc1/Cdc6 (+) \n (+)"),
  c("Orc1/Cdc6 (-)", "Orc1/Cdc6 (+) \n (+)"))

# =======================================
# Step 9: 🎨 Define colors for groups
# =======================================
colors_orcs <- c("Orc1/Cdc6 (-)" = rgb(44, 123, 182, maxColorValue = 255),
                 "Orc1/Cdc6 (+) \n (-) " = rgb(158, 202, 225, maxColorValue = 255),
                 "Orc1/Cdc6 (+) \n (+)"  = rgb(215, 25, 28, maxColorValue = 255))

# =======================================
# Step 10: 🎻 Create plot
# =======================================
p <- ggplot(df_Orcs, aes(x = group, y = SNPsPer1000bp, fill = group)) +
  geom_boxplot(width = 0.6, 
               outlier.shape = NA, 
               alpha = 0.8) +
  geom_jitter(width = 0.1, alpha = 0.3) + # Optional: shows distribution
  theme_classic() +
  labs(y = "SNP density (SNPs/kb)", x= NULL) +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, lineheight = 2),
        axis.text.y = element_text(size = 12),
        legend.position = "none",
        plot.margin = margin(t = 10, r = 10, b = 15, l = 20)) +
  
  scale_fill_manual(values = colors_orcs) +
  geom_signif(comparisons = my_comparisons,
              step_increase = 0.1, # Automatically offsets multiple bars
              map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "ns"=1),
              test = ifelse(is_normal, "t.test", "wilcox.test"), # Individual pair tests
              textsize = 4) +
  ylim(9,15)

q <- p +
  labs(
    y = "SNP Density (SNPs/kb)",
    x = NULL,
    caption = "Fired origins"
  ) +
  theme(
    plot.caption = element_text(hjust = 0.16, vjust = 10 ,  size = 10, margin = margin(t = 10))
  )

# =======================================
# Step 11: 📊 Display & Save
# =======================================
print(q)

output_path_1 <- "~/Documents/DNAscent-Origins/dataAnalisys/Call-Variants/SNPs"
if(dir.exists(output_path_1)) {
  ggsave(filename = file.path(output_path_1, "SNPsInOrc1cdc6Sites_3Groups.pdf"), 
         plot = p, height = 9, width = 18, dpi = 1200, units = "cm")
}

setwd("~/Documents/paper_GenomeVariability/figures/Original/SNPsDensity")
ggsave("SNPsInOrc1cdc6Sites_3Groups.pdf",  plot = q, height = 10.5, width = 16, dpi = 1200, units = "cm")
