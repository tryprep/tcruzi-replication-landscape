# =============================
# Step 1: 📦 Load or install required packages
# =============================
required_packages <- c("ggplot2", "rio", "cowplot", "png", "grid", "viridis")
missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}
suppressMessages(sapply(required_packages, require, character.only = TRUE))

# =============================
# Step 2: 📁 Set working directory
# =============================
setwd("~/Documents/DNAscent-Origins/dataAnalisys/frequency-Origins")

# =============================
# Step 3: 📄 Import and prepare data
# =============================
non_sinc <-  read.csv("mergedPutativeTerminition.bed", quote = "\t", sep = "\t", header = FALSE)
non_sinc$length <- non_sinc$V3 - non_sinc$V2
colnames(non_sinc) <- c("chrom", "start", "end", "length")
df <- rbind(non_sinc)

# =============================
# Step 4: 📊 Compute y_max for histogram
# =============================
binwidth <- 1000
h <- hist(df$length, breaks = seq(0, max(df$length) + binwidth, by = binwidth), plot = FALSE)
y_max <- max(h$counts)

# =============================
# Step 5: 📈 Build plot
# =============================
C <- ggplot(df, aes(x = length)) +
  geom_histogram(binwidth = binwidth, fill = rgb(191,144,0, maxColorValue = 255),
                 color = "black", alpha = 0.8) +
  theme_classic() +
  labs(x = "length of merger (bp)", y = "number of mergers") +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.00000001)
  ) +
  scale_y_continuous(breaks = seq(0, 1400, by = 100), expand = expansion(mult = c(0.01, 0.1))) +
  scale_x_continuous(breaks = seq(0, 20000, by = 2000), expand = expansion(mult = c(0.01, 0.1))) #+
  #annotate("text", x = max(df$length) * 0.3, y = y_max * 0.25, label = "termination", size = 4, fontface = "bold", color = rgb(26, 76, 166, maxColorValue = 255)) +
  #annotate("text", x = max(df$length) * 0.8, y = y_max * 0.25, label = "mini termination zones",  size = 4, fontface = "bold", color = rgb(200, 35, 0, maxColorValue = 255))

print(C)

# =============================
# Step 5: 💾 Save plot
# =============================
# Save in current analysis folder
ggsave("frequencyDnascentTerminationDetermination.pdf", 
       plot = C, height = 9, width = 9, dpi = 1200, units = "cm")

# Save in paper folder
setwd("~/Documents/paper_GenomeVariability/figures/Original/datasetConstrutionandClassification-D-NascentTermination")
ggsave("frequencyDnascentTerminationDetermination.pdf", 
       plot = C, height = 9, width = 18, dpi = 1200, units = "cm")

