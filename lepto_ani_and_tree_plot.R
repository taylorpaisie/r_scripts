rm(list = ls())

# Set the working directory
setwd("/Users/ltj8/Documents/lepto/pr_lepto/redo_all_pr_lepto/ANI_against_type_strain/")

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(viridis)
library(dplyr)
library(ape)
library(ggtree)
library(cowplot)
library(ggplotify)

# Load the mapping between Accessions and Leptospira species names
mapping_df <- read.csv("../lepto_type_strains/lepto_type_strains_ids_and_species.csv")


# Check the results
head(mapping_df)

# Create a list of all "_ANI.Summary.tsv" files
file_list <- list.files(pattern = "_ANI\\.Summary\\.tsv$")

# Initialize an empty dataframe to store combined data
combined_data <- data.frame()

# Loop through the file list, read each file, and combine the data
for (file in file_list) {
  # Read the current file
  temp_data <- read.csv(file, sep = "\t")
  
  # Extract the sample ID from the file name
  sample_id <- gsub("_ANI\\.Summary\\.tsv$", "", file)
  temp_data$SampleID <- sample_id
  
  # Replace 'Bidirectional_ANI...' with the actual column name from your files
  correct_ani_column_name <- 'Bidirectional_ANI...' 
  if (correct_ani_column_name %in% colnames(temp_data)) {
    temp_data <- temp_data[, c('Sample', correct_ani_column_name)]
    colnames(temp_data) <- c('Sample', 'ANI')
    temp_data$SampleID <- sample_id
  } else {
    print(paste("Column", correct_ani_column_name, "not found in file:", file))
    next  # Skip to the next file if the column is not found
  }
  
  # Combine with the main dataframe
  if (nrow(combined_data) == 0) {
    combined_data <- temp_data
  } else {
    combined_data <- rbind(combined_data, temp_data)
  }
}



# Merge the combined data with the mapping dataframe
# Assuming 'Sample' in combined_data corresponds to 'Accession' in mapping_df
combined_data_updated <- merge(combined_data, 
                               mapping_df, by.x = 'Sample', 
                               by.y = 'Accession', all.x = TRUE)

# Check if all samples found a matching species, if not, you might need to address missing matches
missing_species <- combined_data_updated[is.na(combined_data_updated$Leptospira.species), ]
if(nrow(missing_species) > 0) {
  print("Some samples didn't find a matching Leptospira species:")
  print(missing_species$Sample)
}

# Check the column names of combined_data_updated
colnames(combined_data_updated)

if (!"Value" %in% colnames(combined_data_updated)) {
  long_data_updated <- melt(combined_data_updated, 
                            id.vars = c('Leptospira.species', 'SampleID'), 
                            variable.name = 'Metric', 
                            value.name = 'Value')
  
  # Convert 'Value' to numeric to avoid the discrete value error
  long_data_updated$Value <- as.numeric(as.character(long_data_updated$Value))
  
  # Check if there were any NAs introduced by conversion, which indicates non-numeric values were present
  if (any(is.na(long_data_updated$Value))) {
    warning("NAs introduced by conversion indicates non-numeric values in 'Value' column")
  }
  
  # Now, use long_data_updated for the plot
  p_updated <- ggplot(long_data_updated, 
                      aes(x = Leptospira.species, y = SampleID, fill = Value)) +
    geom_tile() +
    # scale_fill_viridis(direction = -1, option = "magma") +
    scale_fill_gradient(low = "#2D00F7", high = "#F20089") +
    labs(title = "ANI Heatmap", x = "Leptospira Species", y = "Sample ID") +
    theme_minimal() +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.title.x = element_text(size=15)) +
    theme(axis.text.y = element_text(size=12)) +
    theme(axis.text.x = element_text(size=16, 
                                     angle = 90, hjust = 1,face = "italic"))
  
  # Display the updated plot
  plot(p_updated)
  
  # Save the updated plot to a file
  ggsave("pr_lepto_ANI_plot_updated.png", plot = p_updated, 
         width = 30, height = 20, bg = "white", dpi = 300)
} else {
  print("The 'Value' column exists. Please check the data structure again.")
}


p_updated <- ggplot(long_data_updated, 
                    aes(x = Leptospira.species, y = SampleID, fill = Value)) +
  geom_tile() +
  # scale_fill_viridis(direction = -1, option = "magma") +
  scale_fill_gradient(low = "#2D00F7", high = "#F20089") +
  labs(title = "ANI Heatmap", x = "Leptospira Species", y = "Sample ID") +
  theme_minimal() +
  theme(axis.title.y = element_text(size=15)) +
  theme(axis.title.x = element_text(size=15)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.text.x = element_text(size=16, 
                                   angle = 90, hjust = 1,face = "italic"))

plot(p_updated)


#####
# starting adding the tree 

tree <- read.tree("../only_pr_lepto/gubbins/gubbins.tree")

# Create a ggtree plot with adjusted margins
tree_plot <- ggtree(tree, 
                    ladderize = TRUE, 
                    layout = "rectangular") +
  geom_tiplab(size=2, color="black", offset = 0.0000005) +
  geom_treescale() #+
# theme(plot.margin = unit(c(0.05,1,1,0.05), "cm"))

plot(tree_plot)

# Convert ggtree plot to ggplot object
tree_plot_gg <- as.ggplot(tree_plot)


tree_plot_gg <- tree_plot_gg +
  theme(plot.margin = margin(t = 25, r = 5, 
                             b = 30, l = 5, unit = "pt"))
ani_plot <- p_updated +
  theme(plot.margin = margin(t = 10, r = 5, 
                             b = 10, l = 5, unit = "pt"))


# Align and combine the tree plot and the recombination plot side by side
combined_plot <- cowplot::plot_grid(tree_plot_gg, 
                                    ani_plot, 
                                    ncol = 2)

plot(combined_plot)


# Save the combined plot to a file
ggsave("combined_tree_ani_plot.png", plot = combined_plot, 
       width = 40, height = 20, bg = "white", dpi = 300)
