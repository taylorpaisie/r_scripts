rm(list = ls())

# setwd('/Users/ltj8/Documents/burk_projects/ms_burk/bayesian_analysis/mol_clock/')
setwd('/Users/ltj8/Documents/burk_projects/2023_october_patient_from_thailand/bayesian_analysis/mol_clock/')

library(ggtree)
library(ape)
library(ggplot2)
library(colorspace)
library(Biostrings)
library(phytools)
library(treeio)
library(dplyr)
library(readr)
library(ggnewscale)

#molecular clock clinical-environmental mcc tree
mcc <- read.beast(file = "rc_bsr/gubbins_recombination_masked_snps_rc_bsr_100_time_mcc.tree")
mlst <- read.csv2(file = "subpanel_mlst_metadata.csv", sep = ',',
                  header = TRUE)
loc <- read.csv2(file = 'subpanel_location_metadata.csv', sep=',',
                 header = TRUE, row.names = 1)
# tree <- read.tree(file = 'OUEST-ONLY-snps_RC_BSP_ASYM_BSVSS_PSSS_MCC.nwk')

get.fields(mcc)

dd <- data.frame(mlst, check.rows = TRUE, check.names = TRUE)
row.names(dd) <- NULL

cb_palette <- c('Novel' = '#5e4fa2ff',
                'NA' = '#1e90bfff',
                '1787' = '#d9434eff',
                '501' = '#bce6a0ff',
                '734' = '#917c6fff',
                '188' = '#75c8a4ff',
                '1034' = '#fdeea1ff',
                '1342' = '#ff7f2aff',
                '1650' = '#008000ff')

loc_colors <- c('Myanmar' = '#fa80ff',
                'Thailand' = '#7eba00',
                'Laos' = '#540054',
                'Malaysia' = '#00dbed',
                'Vietnam' = '#ed3e00',
                'India' = '#004925')



p <- ggtree(mcc, mrsd = "2023-10-01",size=0.75, ladderize=TRUE) +
  geom_point2(aes(label = posterior, subset = posterior >= 0.9),
              size = 3, shape = 18) #+
  # geom_tiplab(size=2, color="black") +
  # geom_tippoint(aes(color=MLST), data = mlst, size = 3.5) +
  # scale_color_manual(values=cb_palette) +
  # theme_tree2()

mcc_tree <- ggtree(mcc, mrsd = "2023-10-01",size=0.75, ladderize=TRUE) +
  # geom_range(range='length_0.95_HPD', color='red', alpha=.6, size=2) +
  geom_point2(aes(label = posterior, subset = posterior >= 0.9),
              size = 3, shape = 18)
plot(mcc_tree)
  
  

# plot(p)
# 
# pp <- p %<+% dd + geom_tippoint(aes(color=MLST), size = 4.5) + 
#   scale_color_manual(values = cb_palette) + 
#   # scale_x_continuous(breaks = seq(2010, 2023)) + 
#   theme_tree2()
# plot(pp)


q <- (p + scale_y_continuous(expand=c(0, 0.3))) %>%
  gheatmap(loc, offset=20, width = 0.05, colnames = FALSE) +
  scale_fill_manual(values = loc_colors, name = "Country") +
  theme_tree2(legend.position='right')
plot(q)

pp <- q %<+% dd + geom_tippoint(aes(color=MLST), size = 4.5) + 
  scale_color_manual(values = cb_palette, name = "MLST") + 
  scale_x_continuous(breaks = c(1825,1850, 1875,1900,1925,1950,1975,2000,2025)) +
  theme_tree2()
plot(pp)

dev.print(pdf, '2023_oct_burk_mcc_tree.pdf')


# 
# 
# p <- ggtree(mcc, aes(color = location), mrsd = "2020-01-01", size = 0.75) + 
#   # scale_color_manual(values=c("mediumpurple4", "springgreen3")) +
#   geom_point2(aes(label = posterior, subset = posterior >= 0.9), size = 2.5) +
#   # scale_x_continuous(breaks = seq(2010, 2023)) + 
#   theme_tree2() + 
#   geom_tiplab(size=2.5, color="black")
# plot(p)







