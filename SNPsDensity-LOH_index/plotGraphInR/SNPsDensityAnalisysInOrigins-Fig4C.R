# ==============================================================================
# SNP DENSITY ANALYSIS PIPELINE - LOCAL SERVER
# Description: This script processes genomic data to analyze SNP density across 
# specific regions (Origins and IZD).
# Input file: SNPs-InOriginsAndTerminationsAndOrcs-Summary.csv
# ==============================================================================

# =======================================
# Step 1: 📦 Load or install required packages
# =======================================
packs <- c("ggplot2", "rio", "dplyr", "plyr", "ggsignif", 
           "rstatix", "RVAideMemoire", "car", "FSA", 
           "dunn.test", "patchwork")

missing_packs <- packs[!packs %in% installed.packages()]
if (length(missing_packs) > 0) {
  install.packages(missing_packs, dependencies = TRUE)
}

sapply(packs, require, character.only = TRUE)

# =======================================
# Step 2: 📂 Set working directory and import dataset
# =======================================
setwd("~/Documents/DNAscent-Origins/dataAnalisys/Call-Variants/SNPs")

data <- read.csv("SNPs-InOriginsAndTerminationsAndOrcs-Summary.csv", 
                 header = TRUE, sep = "\t")

colnames(data) <- c("sample", "count", "length", "group", "SNPsDensity")

# Round SNPsDensity values to 4 decimal places
data$SNPsDensity <- round(data$SNPsDensity, digits = 4)

# =======================================
# Step 3: 🔍 Filter and rename DNAscent-specific groups
# =======================================
DNAscent <- data %>%
  filter(group %in% c("genome-2000-ALL-ORIGIN-FREE",
                      "origins",
                      "miniInitiationZones")) %>%
  mutate(group = dplyr::recode(group,
                               "genome-2000-ALL-ORIGIN-FREE" = "without \n origins",
                               "origins" = "origins",
                               "miniInitiationZones" = "IZD"))

df_DNAscent <- as.data.frame(DNAscent)

# Set factor level order
df_DNAscent$group <- factor(df_DNAscent$group, 
                            levels = c("without \n origins", "origins", "IZD"))
df_DNAscent$SNPsPer1000bp <- df_DNAscent$SNPsDensity * 1000

# =======================================
# Step 4: 📊 Summary statistics
# =======================================
glimpse(df_DNAscent)

mu <- ddply(df_DNAscent, "group", summarise, 
            grp.mean = mean(SNPsDensity, na.rm = TRUE), 
            grp.median = median(SNPsDensity, na.rm = TRUE))

# =======================================
# Step 5: 📈 Test for normality (Shapiro-Wilk)
# =======================================
shapiro_results2 <- sapply(levels(df_DNAscent$group), function(g) {
  group_data <- df_DNAscent$SNPsDensity[df_DNAscent$group == g]
  if (length(group_data[!is.na(group_data)]) > 2) {
    return(shapiro.test(group_data)$p.value)
  } else {
    return(NA)
  }
})

# =======================================
# Step 6: 🧪 Choose and run statistical test
# =======================================
test_type2 <- NULL

if (all(shapiro_results2[!is.na(shapiro_results2)] > 0.05)) {
  anova_result2 <- aov(SNPsDensity ~ group, data = df_DNAscent)
  p_value2 <- summary(anova_result2)[[1]][["Pr(>F)"]][1]
  
  if (p_value2 < 0.05) {
    post_hoc2 <- TukeyHSD(anova_result2)
    test_type2 <- "ANOVA"
  } else {
    signif_labels2 <- rep("ns", 3)
  }
} else {
  kruskal_result2 <- kruskal.test(SNPsDensity ~ group, data = df_DNAscent)
  p_value2 <- kruskal_result2$p.value
  
  if (p_value2 < 0.05) {
    df_DNAscent$SNPsDensity <- jitter(df_DNAscent$SNPsDensity, amount = 1e-5)  # Avoid ties
    post_hoc2 <- dunnTest(SNPsDensity ~ group, data = df_DNAscent, method = "bh")
    test_type2 <- "Kruskal-Wallis"
  } else {
    signif_labels2 <- rep("ns", 3)
  }
}

# =======================================
# Step 7: 🏷️ Prepare comparison annotations
# =======================================
comparison_labels2 <- list(
  c("without \n origins", "origins"),
  c("without \n origins", "IZD"),
  c("origins", "IZD")
)

comparisons_interest2 <- sapply(comparison_labels2, function(x) paste(x, collapse = " - "))

# Extract p-values from post hoc results
if (test_type2 == "Kruskal-Wallis") {
  all_results2 <- post_hoc2$res
  matched_pvals2 <- all_results2[match(comparisons_interest2, all_results2$Comparison), "P.adj"]
} else if (test_type2 == "ANOVA") {
  tukey_df2 <- as.data.frame(post_hoc2$group)
  tukey_df2$Comparison <- rownames(tukey_df2)
  matched_pvals2 <- tukey_df2[match(comparisons_interest2, tukey_df2$Comparison), "p adj"]
}

# Convert p-values into significance symbols
p_to_symbol2 <- function(p) {
  if (is.na(p)) return("ns")
  else if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")
}
signif_labels2 <- sapply(matched_pvals2, p_to_symbol2)

# =======================================
# Step 8: 🎨 Define colors
# =======================================
origin_colors <- c( 
  "without \n origins" = rgb(0, 0, 255,   maxColorValue = 255), 
  "origins"            = rgb(0, 100, 255,  maxColorValue = 255),  
  "IZD"                = rgb(255, 0, 0,    maxColorValue = 255) 
)

# =======================================
# Step 9: 🖼️ Generate Box plot
# =======================================
p <- ggplot(df_DNAscent, aes(x = group, y = SNPsPer1000bp, fill = group)) +
  geom_boxplot(width = 0.8,           # Increase box width (from 0.5 to 0.8)
               linewidth = 0.5,       # Set box line thickness
               outlier.shape = NA, 
               alpha = 1,
               position = position_dodge(width = 0.9)) +
  theme_classic() +
  labs(y = "SNPs density (SNPs/kb)", x = NULL) +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, angle = 0, vjust = 0.7),
        axis.title.y = element_text(vjust = 4, size = 12),
        legend.position = "none",
        plot.margin = margin(t = 10, r = 10, b = 15, l = 20),
        panel.spacing = unit(1, "lines")) +
  scale_fill_manual(values = origin_colors ) +
  geom_signif(comparisons = comparison_labels2,
              annotations = c("ns", "***", "***"),
              y_position = c(13, 14.5, 15.25),
              tip_length = 0.02,
              size = 0.1,        # Increase significance line thickness
              textsize = 5.0)   + 
  ylim(9, 16)

p

# =======================================
# Step 10: 💾 Save plot to desired directories
# =======================================
setwd("~/Documents/DNAscent-Origins/dataAnalisys/Call-Variants/SNPs")
ggsave("SNPsIn-smAtlas.pdf", p, height = 10, width = 10, dpi = 1200, units = "cm")

setwd("~/Documents/paper_GenomeVariability/figures/Original/SNPsDensity")
ggsave("SNPsIn-smAtlas.pdf", p, height = 11, width = 10, dpi = 1200, units = "cm")
