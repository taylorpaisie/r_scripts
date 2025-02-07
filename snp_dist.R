rm(list = ls())

library(reshape2)
library(scales)
library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)

setwd('')

df <- read.csv2('SNP-distances.matrix.tsv', 
                sep = '\t', header = TRUE)
head(df)

# Add a header to the first column
colnames(df)[1] <- "Sample"
head(df)

df_melted <- melt(df, id.vars = "Sample")

# Plot
ggplot(df_melted, aes(x = variable, y = Sample, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c() +
  # scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "SNP Distance Heatmap", x = "Sample", y = "Variable")
