rm(list = ls())

library(reshape2)
library(scales)
library(tidyr)
library(dplyr)
library(ggplot2)

setwd('/Users/ltj8/Documents/lepto/mn_lepto/')

df <- read.csv2('ncbi_lepto_primer_table.csv', 
                sep = ',', header = TRUE)


dfm <- melt(df, value.name = "Count")
head(dfm)

cb_palette <- c('X..of.hits' = '#a670ff',
                'X..of.total.assessed..by.species.' = '#01c47b')



samp_plot <- ggplot(dfm, aes(x = Species,
                            y = Count,
                            fill = variable)) +
  theme_classic() +
  # geom_bar(stat = "identity") +
  geom_bar(position = "dodge",stat = "identity") +
  ggtitle("Leptospira isolates by species") +
  ylab("Number of hits") +  # Change y-axis label
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    legend.title = element_text(vjust = 1, size = 10),
    legend.text = element_text(size = 6),
    axis.title.x = element_text(vjust = 1, size = 16),
    axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
    axis.title.y = element_text(vjust = 1, size = 16),
    axis.text.y = element_text(size = 15)
   ) +
  scale_fill_manual(values = cb_palette, name = "isPcr Hits", 
                    labels = c("Number of hits", "Number of total species assessed"))

plot(samp_plot)

# p <- samp_plot + guides(color = guide_legend(override.aes = list(size = 0.5)))
# plot(p)


dev.print(pdf, 'number_of_ispcr_hits_barplot.pdf')



df <- read.csv2('ncbi_blastn_validation_table.csv', 
                sep = ',', header = TRUE)
head(df)

# Add a header to the first column
# colnames(df)[1] <- "Sample"
# head(df)

dfm <- melt(df, id.vars = "Species")
head(dfm)


cb_palette <- c('Total_._Alignments' = "#01c546",
                'X._Leptospira_Genus_Alignments' = '#00418f',
                'X._Species_Specific_Alignments' = '#ff9b94')



p <- ggplot(dfm, aes(x = Species,
                             y = value,
                             fill = variable)) +
  theme_classic() +
  # geom_bar(stat = "identity") +
  geom_bar(position = "dodge",stat = "identity") +
  ggtitle("Leptospira isolates by species") +
  ylab("Number of hits") +  # Change y-axis label
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    legend.title = element_text(vjust = 1, size = 10),
    legend.text = element_text(size = 6),
    axis.title.x = element_text(vjust = 1, size = 16),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.title.y = element_text(vjust = 1, size = 16),
    axis.text.y = element_text(size = 15)
  ) +
  scale_fill_manual(values = cb_palette, name = "NCBI Leptospira Hits", 
                    labels = c("Number of total alignments",
                               "Number of Leptospira genus alignments",
                               "Number of Leptospira species alignments"))

plot(p)

# p <- samp_plot + guides(color = guide_legend(override.aes = list(size = 0.5)))
# plot(p)


dev.print(pdf, 'number_of_ispcr_hits_barplot.pdf')
