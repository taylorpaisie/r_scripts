rm(list = ls())

setwd("")

library(ggtree)
library(ape)
library(ggsci)
library(ggplot2)
library(colorspace)
library(Biostrings)
library(phytools)
library(treeio)
library(dplyr)
library(readr)
library(phangorn)
library(geiger)
library(RColorBrewer)
library(tangler)

###################################################################
# CA subpanel burk recombination sites masked ML tree

source <- read.csv2(file = "3101811498_CA_subpanel_location.csv", 
                    header = TRUE,
                    sep = ",")
mlst <- read.csv2(file = '3101811498_CA_subpanel_mlst.csv', header = TRUE,
                  sep = ',', row.names = 1)
ml_tree <- read.tree("Parsnp_output/gubbins/gubbins_masked_recombination.fasta.treefile")


q <- ggtree(midpoint(ml_tree), ladderize = TRUE, layout = "rectangular")
plot(q)



d <- q$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label > 90,]


dd <- data.frame(source, check.rows = TRUE, check.names = TRUE)
row.names(dd) <- NULL
dd


cb_palette <- c('USA: CA' = '#d9434eff',
                'Thailand' = '#1e90bfff')#,
                # 'Vietnam' = '#ffbd00')
                # 'USA: MA' = '#bce6a0ff')
# 'USA: IL' = '#ffbd00')
# 'IN' = '#d9434eff',
# 'LK' = '#bce6a0ff',
# 'CN' = '#df5d00',
# 'IL' = '#75c8a4ff',
# 'VN' = '#fdeea1ff')


# # Classic palette BuPu, with 4 colors
# cb_palette <- brewer.pal(5, "Set1")
# 
# # Add more colors to this palette :
# cb_palette <- colorRampPalette(cb_palette)(4)

p <- q %<+% source + 
  geom_nodepoint(data=d, aes(label=label), size = 2.5, shape = 18) +
  geom_tippoint(aes(color=Location), size = 3.75) +
  geom_tiplab(size=3, color="black", offset = 0.00000002) +
  scale_color_manual(values = cb_palette) +
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_startrek() +
  # scale_color_viridis(discrete = TRUE) +
  geom_treescale(y=-2, offset=0.5) +
  theme(legend.position = c(0, 1), 
        legend.justification = c(0, 1))
p

mlst_color <- c('ST-708' = '#ef476f',
  'ST-212' = '#26547c',
  'ST-211' = '#ffd166')
# 
# 
pp <- (p + scale_y_continuous(expand = c(0, 0.2))) %>%
   gheatmap(mlst, offset = 0.001, width = 0.05,
            legend_title = "MLST", colnames = FALSE) +
   scale_fill_manual(values = mlst_color, name = "MLST")

pp

# # Create the plot with Phylopic images
# p <- q %<+% dd +
#   geom_tippoint(aes(color=Species), size = 3) +
#   geom_nodepoint(data=d, aes(label=label), size = 2.5, shape = 18) +
#   geom_tiplab(aes(image=image_uuid),
#               geom="phylopic") +
#   # geom_image(aes(image=image_url), size=0.05, offset=1) +
#   scale_color_manual(values = cb_palette) +
#   theme(legend.position = c(0.1, 1.0),
#         legend.justification = c(0, 1)) +
#   geom_treescale()
#
# p


# dev.print(pdf, 'ma_burk_subpanel_ml_tree.pdf')
dev.print(pdf, '3101811498_CA_subpanel_mlst_ml_tree.pdf', 
          width = 13, height = 8.5)

ggsave("3101811498_CA_subpanel_mlst_ml_tree.png", plot = p, 
       # width = 20, height = 20, 
       dpi = 300)
