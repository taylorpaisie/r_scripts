rm(list = ls())

setwd('/Users/ltj8/Documents/burk_projects/ms_burk/subpanel/snps')

# data.folder <- file.path('/Users/ltj8/Documents/burk_projects/2023_october_patient_from_thailand/bayesian_analysis/snps_iqtree/')

library(ggtree)
library(ape)
library(ggplot2)
library(colorspace)
library(Biostrings)
library(phytools)
library(treeio)
library(dplyr)
library(readr)
library(phangorn)
library(geiger)


source <- read.csv2(file = "ms_burk_subpanel_source_metadata.csv", header = TRUE, 
                    sep = ",")
# mlst <- read.csv2(file = 'ms_burk_mlst_metadata.csv', header = TRUE,
#                   sep = ',', row.names = 1)
ml_tree <- read.tree("../gubbins/gubbins.tree")


q <- ggtree(midpoint(ml_tree), ladderize = TRUE, layout = "rectangular")
plot(q)

d <- q$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label > 90,]


dd <- data.frame(source, check.rows = TRUE, check.names = TRUE)
row.names(dd) <- NULL

cb_palette <- c('Clinical' = 'mediumpurple4', 
                'Environment' = 'springgreen3')
# cb_palette <- c('Novel' = '#5e4fa2ff',
#                 'NA' = '#1e90bfff',
#                 'ST-1787' = '#d9434eff',
#                 'ST-501' = '#bce6a0ff',
#                 'ST-734' = '#917c6fff',
#                 'ST-188' = '#75c8a4ff',
#                 'ST-1034' = '#fdeea1ff',
#                 'ST-1342' = '#ff7f2aff',
#                 'ST-1650' = '#008000ff')

p <- q %<+% dd + geom_tippoint(aes(color=Source), size = 4.5) +
  geom_nodepoint(data=d, aes(label=label), size = 2.5, shape = 19) + 
  geom_tiplab(size=2, color="black", offset = 0.00005) +
  scale_color_manual(values = cb_palette) +
  geom_treescale()

  # geom_strip(54, 66, label = '2014-2015', align = T, color = 'pink', barsize = 3) +
  # geom_strip(83, 53, align = T, label = '2013', color = '#3FA9F5', barsize = 3) +
  # geom_strip(111, 91, label = '2010-2012', align = T, color = 'darkgrey', barsize = 3) + 
  # geom_strip(34, 87, label = '2012-2013', align = T, color = 'gold', barsize = 3) + 
  # geom_hilight(node = 123, fill = "pink", alpha = 0.2) +
  # geom_hilight(node = 169, fill = 'darkgray', alpha = 0.2) + 
  # geom_hilight(node = 150, fill = '#3FA9F5', alpha = 0.2) + 
  # geom_hilight(node = 204, fill = 'gold', alpha = 0.2) + 
plot(p)

mlst_color <- c("ST-92" = "#1e90bfff",
                "ST-544" = '#ff7f2aff')

# country_color <- c("Myanmar" = "", "Thailand", "Malaysia","Vietnam", "India")  

pp <- (p + scale_y_continuous(expand=c(0, 0.3))) %>% 
  gheatmap(mlst, offset=0, width = 0.05, colnames = FALSE) + 
  scale_fill_manual(values = mlst_color, name = "MLST") +
  theme_tree2(legend.position='right')

plot(pp)

dev.print(pdf, 'ms_burk_ml_tree.pdf')
# dev.print(png, 'ms_burk_ml_tree.png')





