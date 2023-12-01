rm(list = ls())

detachAllPackages <- function() {
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}

detachAllPackages()

install_load <- function(Required_Packages) {
  for(package in Required_Packages){
    if (!package %in% installed.packages()) install.packages(package, character.only = TRUE)
    library(package, character.only = TRUE)
  }
}

library(gplots)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(grid)
library(lattice)
library(gridExtra)
library(ComplexHeatmap)

setwd('/Users/ltj8/Documents/lepto/pr_lepto/redo_all_pr_lepto/')


table_input <- read.csv("borgpetersenii/roary_results/gene_presence_absence.csv", 
                        sep=",", na.strings=c("","NA"))
table_input <- as.data.frame(table_input)

table_values <- within(table_input, rm("Gene",
                                       "Non.unique.Gene.name","No..isolates",
                                       "No..sequences",
                                       "Avg.sequences.per.isolate",
                                       "Genome.Fragment","Order.within.Fragment",
                                       "Accessory.Fragment","Accessory.Order.with.Fragment",
                                       "QC","Min.group.size.nuc","Max.group.size.nuc",
                                       "Avg.group.size.nuc"))

abscence_presence <- as.matrix(table_values[,-1])
rownames(abscence_presence) <- table_values[,1]
abscence_presence[is.na(abscence_presence)] <- 0
abscence_presence[which(abscence_presence!=0)] <- 1


a_p_matrix <- mapply(abscence_presence, FUN=as.numeric)
a_p_matrix <- matrix(data=a_p_matrix, ncol=length(colnames(abscence_presence)), 
                     nrow=length(row.names(abscence_presence)))
row.names(a_p_matrix) <- row.names(abscence_presence)
colnames(a_p_matrix) <- colnames(abscence_presence)

my_colors <- c("gray", "#01c47b")

# 
Heatmap(
  a_p_matrix,
  heatmap_legend_param = list(
    title = "Absence/Presence of Genes",
    at = c(0,1),
    labels = c("Absent", "Present")
  ),
  name = "Absence/Presence of Genes",
  col = my_colors,
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 10)
)



# heatmap.2(a_p_matrix,
#           col = c("red", "darkblue"),
#           main = "Absence/Presence of genes",
#           trace = "none",
#           labRow = FALSE,
#           cex.axis = 0.8)  # Adjust the value (0.8 in this case) for the desired font size


