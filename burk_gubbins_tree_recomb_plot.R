rm(list = ls())

library(gplots)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(grid)
library(lattice)
library(gridExtra)
library(ComplexHeatmap)
library(ape)
library(ggtree)
library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(viridis)
library(cowplot)
library(ggplotify)

setwd('/Users/ltj8/Documents/burk_projects/ms_burk/new_subpanel/')


# Load the GFF file
gff_file_path <- "gubbins/gubbins.recombination_predictions.gff"

# Read the GFF file
gff_data <- read.delim(gff_file_path, 
                       header = FALSE, comment.char = '#')

# Select only the relevant columns and rename them
recombination_data <- gff_data %>%
  select(Taxa = V9, Start = V4, End = V5) %>%
  mutate(Taxa = sub('.*;taxa=([^;]+).*', '\\1', Taxa),
         Start = as.numeric(Start),
         End = as.numeric(End))

# List of taxa to include
taxa_list <- c("3001284746", "3001964113", "3002023864", "3002023865", "3002023871",
    "3002023874", "3002023876", "3002023878", "3003781053", "3003781054",
    "3003787861", "3003787862", "3003787863", "3003787894", "3003787895",
    "3003787897", "3003787898", "3003787899", "3003787907", "3003787912",
    "3003787937", "GCF_000959265.1", "GCF_001885195.1", "GCF_002110925.1",
    "GCF_002110945.1", "GCF_002110965.1", "GCF_002110985.1", "GCF_002111005.1",
    "GCF_002111025.1", "GCF_002111045.1", "GCF_002111065.1", "GCF_002111085.1",
    "GCF_002111105.1", "GCF_002111125.1", "GCF_002111145.1", "GCF_002111165.1",
    "GCF_002111185.1", "GCF_002111205.1", "GCF_002111245.1", "GCF_002111265.1",
    "GCF_002111285.1", "GCF_002111305.1", "GCF_002111325.1", "GCF_002111345.1",
    "GCF_002111385.1", "GCF_002113945.1", "GCF_002115385.1", "GCF_003583425.1",
    "GCF_003583435.1", "GCF_003584055.1", "GCF_003584065.1", "GCF_003813825.1",
    "GCF_003813945.1", "GCF_003854325.1", "GCF_003858035.1", "GCF_003944665.1",
    "GCF_006542565.1", "GCF_006542585.1", "GCF_007995115.1", "GCF_013265625.1",
    "GCF_013265665.1", "GCF_013265695.1", "GCF_015714675.1", "GCF_016092155.1",
    "GCF_016617505.1", "GCF_023383335.1")

tree_file <- "gubbins/gubbins.tree"
phylo_tree <- read.tree(tree_file)

# Assuming phylo_tree is your phylogenetic tree object
# Extract the order of taxa from the phylogenetic tree
tree_taxa <- phylo_tree$tip.label

# Reorder taxa_list to match the order of taxa in the phylogenetic tree
ordered_taxa_list <- intersect(tree_taxa, taxa_list)

# Ensure the recombination data is filtered and ordered according to the tree
filtered_recombination_df <- recombination_data %>%
  filter(Taxa %in% ordered_taxa_list) %>%
  mutate(Taxa = factor(Taxa, levels = ordered_taxa_list))

# # Filter the data for the taxa list
# filtered_recombination_df <- recombination_data %>%
#   filter(Taxa %in% taxa_list)

# Create a mapping of taxa to y-axis values
taxa_to_y <- setNames(seq_along(taxa_list), taxa_list)

# Creating the plot
recomb_plot <- ggplot() + 
  geom_rect(data = filtered_recombination_df, 
            aes(xmin = Start, xmax = End, 
                ymin = Taxa, ymax = Taxa,
                fill = Taxa), 
            color = "black", size = 3, alpha = 0.75) +
  scale_y_discrete(limits = taxa_list) +
  scale_fill_viridis_d() +
  labs(x = 'Sequence Position', y = 'Taxa',
       title = 'Recombination Events for Provided Taxa') +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(size=15)) +
  theme(axis.title.x = element_text(size=15)) +
  theme(axis.text.y = element_text(size=7)) +
  theme(axis.text.x = element_text(size=10))

plot(recomb_plot)

# Save the plot
# ggsave('recombination_plot_enhanced.png', width = 10, height = 6, dpi = 300)

# Create a ggtree plot with adjusted margins
tree_plot <- ggtree(phylo_tree, 
                    ladderize = TRUE, 
                    layout = "rectangular") +
  geom_tiplab(size=2, color="black", offset = 0.0000005) +
  geom_treescale() #+
# theme(plot.margin = unit(c(0.05,1,1,0.05), "cm"))

# Convert ggtree plot to ggplot object
tree_plot_gg <- as.ggplot(tree_plot)


tree_plot_gg <- tree_plot_gg +
  theme(plot.margin = margin(t = 25, r = 5, 
                             b = 30, l = 5, unit = "pt"))
recomb_plot <- recomb_plot +
  theme(plot.margin = margin(t = 10, r = 5, 
                             b = 10, l = 5, unit = "pt"))


# Align and combine the tree plot and the recombination plot side by side
combined_plot <- cowplot::plot_grid(tree_plot_gg, 
                                    recomb_plot, 
                                    ncol = 2)

plot(combined_plot)

# Save the combined plot to a file
ggsave("combined_tree_recombination_plot.png", plot = combined_plot, 
       width = 40, height = 20, bg = "white", dpi = 300)




