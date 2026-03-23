# This script creates bar plots to visualize DNA replication origin data
# using ggplot2, including a Kruskal-Wallis test.

# First, install and load the required packages
packs <- c("ggplot2", "rio", "cowplot" , "png", "grid", "viridis")

# Check and install missing packages
if (sum(as.numeric(!packs %in% installed.packages())) != 0) {
  installer <- packs[!packs %in% installed.packages()]
  for (i in 1:length(installer)) {
    install.packages(installer, dependencies = TRUE)
    break()
  }
  sapply(packs, require, character.only = TRUE)
} else {
  sapply(packs, require, character.only = TRUE)
}

# Adjust the working directory
setwd("~/Documents/DNAscent-Origins/dataAnalisys/Length")

# Insert the data for the data frame composition
NonSinc   <- read.csv2('StatsLengthRepliconSinc.csv')

# Renaming to columns
colnames(NonSinc)  <- c('Length','Count')

# dataframe composition (database)
df<-as.data.frame(NonSinc)
df
df$Length <- factor(df$Length, levels=c("1-500", "501-1000", "1001-1500", "1501-2000", "2001-2500",
                                        "2501-3000","3001-3500","3501-4000","4001-4500",
                                        "4501-5000","5001-6000", "6001-7000", "7001-8000", "8001-9000", 
                                        "9001-10000" ,"≥ 10.001"))
df
# plot graph

B <- ggplot(df, aes(x = Length, y = Count, fill = as.factor(Length))) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d(option = "turbo") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = -0.01), 
        axis.title.y = element_text(size = 10, angle = 90, vjust = -0.00001),
        axis.title.x = element_text(size = 10, vjust = -0.00001), 
        axis.text = element_text(size = 10, colour = "Black")) +
  coord_cartesian(ylim = c(0, 2000)) +
  scale_y_continuous(breaks = seq(0, 2000, by = 200), expand = expansion(mult = c(0.005, 0.05))) +
  theme(legend.position = "none") +
  labs(x = "Inter-fork distance (bp)", y = "Count", fill = "Length") 
B

ggsave("LenghtPutativeOrigins.pdf", B, height=9, width=9, dpi=1200, , units = "cm")

# set directories to figure 
setwd("~/Documents/paper_GenomeVariability/figures/Original/datasetConstrutionandClassification-D-NascentOrigins")
ggsave("LenghtPutativeOrigins.pdf", B, height=9, width=18, dpi=1200, , units = "cm")

