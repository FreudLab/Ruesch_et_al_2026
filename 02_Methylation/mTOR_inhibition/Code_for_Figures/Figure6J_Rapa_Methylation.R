#Load libraries
library(RnBeads)
library(RnBeads.hg38)
library(ggrepel)
library(tidyverse)
library(GO.db)

#Load the RnBeads Set for the Rapamycin Study (Illumina Epic v2 for this study)
#This RnBeads set also included some other data we are labeling AML
Rapa_RnB <- load.rnb.set(
  "~/DNAmethylation/RnBeads_Output/AML_Rapa_RnBeads/reports/rnbSet_preprocessed"
  )
rapa_pheno <- pheno(Rapa_RnB)
rapa_pheno$CellType <- c(rep("AML", 14), 
                         rep(c("CD16Neg_15D", "CD16Neg_DMSO", "CD16Neg_RAPA"),3),
                         rep(c("CD16Pos_15D", "CD16Pos_DMSO","CD16Pos_RAPA"),3)
)
slot(Rapa_RnB, "pheno") <- rapa_pheno

#-----------------------------------------------------------------------------
#Loop to create each diffmeth object
# 1. Define all your comparisons in a single data frame for clarity
comparisons <- data.frame(
  group1 = c("CD16Neg_15D", "CD16Neg_RAPA", "CD16Pos_15D", "CD16Pos_RAPA"),
  group2 = c("CD16Neg_DMSO", "CD16Neg_DMSO", "CD16Pos_DMSO", "CD16Pos_DMSO"),
  stringsAsFactors = FALSE
)

# 2. Define the base path for saving your files
output_dir <- "~/DNAmethylation/NKDI_NK_IV_Analysis/Rapamycin/"

# 3. Loop through each row of the comparisons data frame
for (i in 1:nrow(comparisons)) {
  
  # Get the two cell types for the current comparison
  type1 <- comparisons$group1[i]
  type2 <- comparisons$group2[i]
  
  # Create a descriptive name for the comparison (e.g., "CD16Neg_15DvDMSO")
  comp_name <- paste0(sub("_15D|_RAPA", "", type1), "_", sub(".*_", "", type1), "v", sub(".*_", "", type2))
  
  # --- Start of your original logic, now using variables ---
  
  cat("Starting comparison:", comp_name, "\n") # Print progress
  
  # Find indices of samples to REMOVE (the ones not in our current comparison)
  remove_idx <- which(!pheno(Rapa_RnB)$CellType %in% c(type1, type2))
  
  # Create the subsetted RnB object by removing those samples
  rnb_subset <- remove.samples(Rapa_RnB, remove_idx)
  
  # Check and print the samples in the current subset
  cat("Samples in subset:\n")
  print(table(pheno(rnb_subset)$CellType))
  
  # Run the differential methylation analysis
  diff_meth_result <- rnb.execute.computeDiffMeth(
    rnb_subset,
    pheno.cols = "CellType"
  )
  
  # Construct the full file path for saving the result
  output_filename <- file.path(output_dir, paste0(comp_name, "_diffmeth"))
  
  # Save the result
  save.rnb.diffmeth(diff_meth_result, output_filename)
  
  cat("Successfully saved results to:", output_filename, "\n\n")
}

#----------------------------------------------------------------------------
#Perform GO Enrichment
CD16Neg_15D_GO <- performGoEnrichment.diffMeth(Rapa_RnB, CD16Neg_15D)
CD16Pos_15D_GO <- performGoEnrichment.diffMeth(Rapa_RnB, CD16Pos_15D)
CD16Neg_RAPA_GO <- performGoEnrichment.diffMeth(Rapa_RnB, CD16Neg_Rapa)
CD16Pos_RAPA_GO <- performGoEnrichment.diffMeth(Rapa_RnB, CD16Pos_Rapa)
#If saving GO, save it as an RDS file for loading later

#Graphing all GO data results

# Create a named list of the key results to iterate over
# We'll look at the top 100 hypomethylated genes in Biological Process (BP)
go_results_list <- list(
  "CD16Neg_15D" = CD16Neg_15D_GO$region[[1]]$BP$genes$rankCut_500$hypo,
  "CD16Pos_15D" = CD16Pos_15D_GO$region[[1]]$BP$genes$rankCut_500$hypo,
  "CD16Neg_RAPA" = CD16Neg_RAPA_GO$region[[1]]$BP$genes$rankCut_500$hyper,
  "CD16Pos_RAPA" = CD16Pos_RAPA_GO$region[[1]]$BP$genes$rankCut_500$hyper
)

# Use purrr::map_dfr to loop through the list, summarize, and combine
# .id = "comparison" creates a new column with the names from our list
all_go_results <- map_dfr(go_results_list, ~summary(.), .id = "comparison")

# Inspect the combined data frame
head(all_go_results)

# Calculate -log10(Pvalue) for plotting
all_go_results <- all_go_results %>%
  mutate(log10Pval = -log10(Pvalue))

# Find the top 7 GO terms from each comparison group
top_terms <- all_go_results %>%
  group_by(comparison) %>%
  slice_max(order_by = log10Pval, n = 5) %>%
  ungroup()

# Get a unique list of these top terms to plot
unique_top_terms <- unique(top_terms$GOBPID)

# Filter the full results to include only these top terms
plot_data <- all_go_results %>%
  filter(GOBPID %in% unique_top_terms)

plot_data$GO_Term <- Term(GOTERM[plot_data$GOBPID])

p1 <- ggplot(plot_data, aes(x = comparison, y = fct_reorder(GO_Term, log10Pval))) +
  geom_point(aes(color = log10Pval, size = Count)) +
  scale_color_gradient(low = "lightgrey", high = "darkred") +
  scale_size(range = c(2, 8)) + # Adjust the range for dot sizes
  theme_bw() + # A clean theme
  labs(
    title = "GO Enrichment for Hypomethylated Regions",
    subtitle = "Top Biological Processes",
    x = "Comparison to 15-C DMSO",
    y = "Gene Ontology Term",
    color = "-log10(P-value)",
    size = "Gene Count"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11), # Rotate x-axis labels
    axis.text.y = element_text(size = 11) # Adjust y-axis font size if needed
  )

p1