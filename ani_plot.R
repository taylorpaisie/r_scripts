rm(list = ls())

# Set the working directory
setwd("")

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(viridis)
library(dplyr)
library(ape)
library(ggtree)
library(gplots)
library(cowplot)
library(ggplotify)
library(dendextend)
library(gplots)
library(RColorBrewer)

par(cex.main=0.72)

# color palette and breaks
# custom_pallette <- colorRampPalette(c('blue', 'red'))
myCol <- c('grey99', 'grey88', 'gray50', 'black')
## Defining breaks for the color scale
myBreaks <- c(90, 95, 97, 99, 100)

# read in data
ani=read.table('ANI.matrix.twoway.tab', header=T, 
               sep='\t', row.names=1, check.names=FALSE)

# Classic palette BuPu, with 4 colors
coul <- brewer.pal(8, "Set3") 

# Add more colors to this palette :
coul <- colorRampPalette(coul)(8)


ani.dist=dist(as.matrix(ani), method ='euclidean')
ani.hclust=hclust(ani.dist, method ='ward.D2')
ani.dnd <- as.dendrogram(ani.hclust)
# cols_branches <- c('black', 'red', 'red4', 'darkgreen', 'dodgerblue3', 'navy')
ani.dnd <- color_branches(ani.dnd, k=8, col=coul)
col_labels <- get_leaves_branches_col(ani.dnd)
col_labels <- col_labels[order(order.dendrogram(ani.dnd))]

# # Remove everything after and including "_GCA" in row names
# rownames(ani) <- make.unique(gsub("_GCA.*", "", rownames(ani)))
# 
# # Remove everything after and including "_GCA" in column names
# colnames(ani) <- make.unique(gsub("_GCA.*", "", colnames(ani)))

# Remove everything after and including the first underscore, but only if there's no "_GCA" in the name
# rownames(ani) <- ifelse(grepl("_GCA", rownames(ani)),
#                         rownames(ani),
#                         gsub("_.*", "", rownames(ani)))
# 
# colnames(ani) <- ifelse(grepl("_GCA", colnames(ani)),
#                         colnames(ani),
#                         gsub("_.*", "", colnames(ani)))



# rownames(ani) <- gsub("_", " ", rownames(ani))
# colnames(ani) <- gsub("_", " ", colnames(ani))

# print sample names to a csv file
# Write matrix to CSV with row names
# write.csv(ani, file = "matrix_output.csv", row.names = TRUE)


# Random values for each Bacillus type
Values <- runif(27, 50, 100)
# 
# Names of Bacillus types, modified for superscript
# Bacillus <- c('Bacillus pacificus EB422 T',
#               'Bacillus pseudomycoides DSM 12442 T',
#               'Bacillus cytotoxicus NVH 39110098 T',
#               'Bacillus thuringiensis ATCC 10792 T',
#               'Bacillus bingmayongensis FJAT10013831 T',
#               'Bacillus arachidis SY8 T',
#               'Bacillus mobilis 0711P91001 T',
#               'Bacillus rhizoplanae CIP 111899 T',
#               'Bacillus anthracis ATCC 14578 T',
#               'Bacillus mycoides DSM 2048 T',
#               'Bacillus sanguinis BML100BC004 T',
#               'Bacillus toyonensis BCT1007112 T',
#               'Bacillus albus N35100101002 T',
#               'Bacillus paramycoides NH24A2 T',
#               'Bacillus weihenstephanensis DSM 11821 T',
#               'Bacillus proteolyticus TD42 T',
#               'Bacillus nitratireducens 4049 T',
#               'Bacillus tropicus N24 T',
#               'Bacillus hominis BML100BC059 T',
#               'Bacillus paranthracis MN5 T',
#               'Bacillus cereus ATCC 14579 T',
#               'Bacillus wiedmannii FSL W81000169 T',
#               'Bacillus manliponensis BL41006 T',
#               'Bacillus paramobilis BML100BC017 T',
#               '3101410416',
#               'Bacillus luti TD41 T',
#               'Bacillus gaemokensis BL31006 T')

# Convert names to expressions for superscript formatting
# Bacillus_labels <- sapply(Bacillus, function(x) {
#   if (grepl(" T$", x)) {
#     # If ' T' is at the end of the name, format it as a superscript
#     main_name <- sub(" T$", " ", x)  # Remove the ' T' from the end
#     # Create the correct expression with superscript
#     as.expression(substitute(paste(XXX, " "^T), list(XXX = main_name)))
#   } else {
#     # If no ' T' at the end, return the name as is
#     x
#   }
# })

# Bacillus_labels <- sapply(Bacillus, function(x) {
#   if (grepl(" T$", x)) {
#     # If ' T' is at the end of the name, format it as a superscript
#     main_name <- sub("T$", "", x)  # Remove the ' T' from the end
#     # Create the correct expression with superscript
#     as.expression(substitute(paste(XXX, ""^T), list(XXX = main_name)))
#   } else {
#     # If no ' T' at the end, return the name as is
#     x
#   }
# })
# # Print to verify the correct formatting
# print(Bacillus_labels[1:5])

# Assuming Values and Bacillus_labels are already defined as before
# Displaying the plot
# barplot(Values, names.arg = Bacillus_labels, las=2, cex.names=0.7, main="Bacillus Values", ylim=c(0, max(Values) + 10))


# Apply the function to column names
# colnames(ani) <- superscript_T(colnames_ani)
# new_colnames

# colnames(ani) <- Bacillus_labels



# heatmap.2(as.matrix(ani), col=myCol, distfun=dist,
#           density.info=c('none'), key.title='',
#           key.xlab='Identity [%]', key.ylab='',
#           breaks=myBreaks, trace='none', margins=c(11, 12),
#           cexCol=0.2 + 1/log10(2 * ncol(ani)),
#           cexRow=0.2 + 1/log10(2 * nrow(ani)),
#           main='Average Nucleotide Identity', key=FALSE,
#           RowSideColors=col_labels, colRow=col_labels,
#           srtCol=45, labCol = Bacillus_labels, labRow = Bacillus_labels)


heatmap.2(as.matrix(ani), col=myCol, distfun=dist,
          density.info=c('none'), key.title='',
          key.xlab='Identity [%]', key.ylab='',
          breaks=myBreaks, trace='none', margins=c(11, 12),
          cexCol=0.2 + 1/log10(2 * ncol(ani)),
          cexRow=0.2 + 1/log10(2 * nrow(ani)),
          main='Average Nucleotide Identity', key=FALSE,
          RowSideColors=col_labels, colRow=col_labels,
          srtCol=45)

legend(x=0.0001, y=1.05, xjust=0.1, cex=0.7, fill=myCol,
       title='Bidirectional Identity [%]',
       inset=c(0.0005,0.005),
       legend=c('< 95', '95 - 97', '98 - 99', '> 99'))

# heatmap.2(as.matrix(ani), col=myCol, distfun=dist,
#           density.info=c('none'), key.title='',
#           key.xlab='Identity [%]', key.ylab='',
#           breaks=myBreaks, trace='none', margins=c(11, 12),
#           cexCol=0.05 + 1/log10(2 * ncol(ani)),
#           cexRow=0.05 + 1/log10(2 * nrow(ani)),
#           main='Average Nucleotide Identity', key=FALSE,
#           RowSideColors=col_labels, colRow=col_labels,
#           srtCol=45)
# 
# 
# legend(x=-0.03, y=1.12, xjust=0, cex=0.65, fill=myCol,
#        title='Bidirectional Identity [%]',
#        inset=c(0.0005,0.005),
#        legend=c('< 95', '95 - 97', '98 - 99', '> 99'))



dev.print(pdf, 'pomona_only_ani_plot.pdf', 
          width = 12, height = 8.5)


