rm(list = ls())

library(reshape2)
library(scales)
library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)

setwd('/Users/ltj8/Documents/lepto/pr_lepto/redo_all_pr_lepto/')

df <- read.csv2('interrogans/parsnp/SNP-distances.matrix.tsv', 
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

df2 <- read.csv2('borgpetersenii/parsnp/SNP-distances.matrix.tsv', 
                sep = '\t', header = TRUE)
head(df2)

# Add a header to the first column
colnames(df2)[1] <- "Sample"
head(df2)

df2_melted <- melt(df2, id.vars = "Sample")

# Plot
ggplot(df2_melted, aes(x = variable, y = Sample, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c() +
  # scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "SNP Distance Heatmap", x = "Sample", y = "Variable")

df3 <- read.csv2('kirscheri/parsnp/SNP-distances.matrix.tsv', 
                 sep = '\t', header = TRUE)
head(df3)

# Add a header to the first column
colnames(df3)[1] <- "Sample"
head(df3)

df3_melted <- melt(df3, id.vars = "Sample")

# Plot
ggplot(df3_melted, aes(x = variable, y = Sample, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "inferno") +
  # scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "SNP Distance Heatmap", x = "Sample", y = "Variable")
