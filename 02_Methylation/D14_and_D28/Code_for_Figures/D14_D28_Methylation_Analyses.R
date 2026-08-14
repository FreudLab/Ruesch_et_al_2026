#Load Packages
library(RnBeads)
library(RnBeads.hg38)
library(ggrepel)
library(tidyverse)

#----------------------------------------------------------------------------
#Necessary for all downstream analyses!
#Load in RnBeads Set and Confirm Phenotype of Samples
#Add in CellType value in phenotype to help with analysis
NK_RnB <- load.rnb.set(
  "~/DNAmethylation/RnBeads_Output/NKDI_and_IV_RnBeadsSet/reports/rnbSet_preprocessed"
  )
phenotype <- pheno(NK_RnB)
phenotype$CellType <- c(rep("15D_D14", 2), 
                        rep("15C_D14", 2), 
                        rep(c("15D_CD16Neg", "15D_CD16Pos"),2),
                        rep(c("15C_CD16Neg", "15C_CD16Pos"), 2),
                        rep("15C_D14", 2),
                        rep(c("15D_CD16Neg", "15D_CD16Pos"),2),
                        "15C_CD16Neg", "15C_CD16Pos", "15C_CD16Neg",
                        rep("S1",4), 
                        rep("S2A", 4),
                        rep("S2B", 4),
                        rep("S3",4), 
                        rep("S4A", 4),
                        rep("S4B", 4),
                        rep("S5", 3)
)
slot(NK_RnB, "pheno") <- phenotype

#----------------------------------------------------------------------------
###Figure 2F and Figure S5D Code###
###PCA of NK Developmental Score based of the eNKTL paper method looking at###
###Stage 3 Tonsil to Stage 6 Blood developmentally variable probes###
# Extract Beta values from NK_RnB
meth_data <- meth(NK_RnB, type = "sites", row.names = TRUE)
colnames(meth_data) <- phenotype$Sample_ID

#Read in 805K preferred probes dataframe and subset probeset on those probes
pref_probes <- read.csv("~/DNAmethylation/NKDI_NK_IV_Analysis/805KAnnotation.csv")
filtered_meth <- meth_data[rownames(meth_data) %in% pref_probes$IlmnID, ] 
# Ensure row names in meth_data correspond to IlmnID
filtered_meth <- filtered_meth[complete.cases(filtered_meth), ]

#Identify ex vivo NKDIs from the data and subset on those samples
BT_RnB <- load.rnb.set(
  "~/DNAmethylation/RnBeads_Output/B_T_IV/reports/rnbSet_preprocessed"
  )
meth_BT <- meth(BT_RnB, type = "sites", row.names = TRUE)
phenotype_BT <- pheno(BT_RnB)
filtered_BT <- meth_BT[rownames(meth_BT) %in% pref_probes$IlmnID, ] 
# Ensure row names in meth_BT correspond to IlmnID
filtered_BT <- filtered_BT[complete.cases(filtered_BT), ]
ex_vivo_meth <- filtered_BT[, 34:54]

#Find the top 5000 variable probes
probe_variances <- apply(ex_vivo_meth, 1, var, na.rm = TRUE)
top5000_probes <- names(sort(probe_variances, decreasing = TRUE)[1:5011])
#Used 5011 so that when filtering on valid probes it equals 5000

#Subset the full dataset on the 5000 probes
valid_probes <- top5000_probes[top5000_probes %in% rownames(filtered_meth)]
pca_data <- na.omit(filtered_meth[valid_probes, ])

#Perform PCA and plot
pca_result <- prcomp(t(pca_data), scale. = TRUE)

pca_df <- data.frame(pca_result$x, 
                     CellType = phenotype$CellType, 
                     SampleID = rownames(phenotype)
                     )
#Save pca_df as csv for importing data and figure generation in graphpad
rownames(pca_df) <- phenotype$Sample_ID
write.csv(pca_df, file='~/DNAmethylation/NKDI_NK_IV_Analysis/PCA_Values.csv')
#----------------------------------------------------------------------------
###Figure 2G and Figure S5C Code###
###epiCMIT scoring of proliferation among D28 cells###
library(GenomicRanges)
library(data.table)

#Load in .RData file for performing epiCMIT proliferation estimate
#Instructions at https://duran-ferrerm.github.io/Pan-B-cell-methylome/Estimate.epiCMIT.html#3_Analyses
load("~/epiCMIT/Estimate.epiCMIT.RData")

beta_vals <- meth(NK_RnB, type='sites', row.names = T)
colnames(beta_vals) <- phenotype$Sample_ID
beta_vals <- as.data.frame(beta_vals)

DNAm.epiCMIT <- DNAm.to.epiCMIT(DNAm = beta_vals,
                                DNAm.genome.assembly = "hg38",
                                map.DNAm.to = "Illumina.450K.epiCMIT",
                                min.epiCMIT.CpGs = 800 # minimum recommended
)

##calculate epiCMIT
epiCMIT.Illumina <- epiCMIT(DNAm.epiCMIT = DNAm.epiCMIT,
                            return.epiCMIT.annot = FALSE,
                            export.results = FALSE#,
                            #export.results.dir = ".",
                            #export.results.name = "Illumina.450k.example_"
)

#Create csv file 
write.csv(epiCMIT.Illumina$epiCMIT.scores, 
          "~/DNAmethylation/NKDI_NK_IV_Analysis/R_Files/epicmit_scores.csv")
#--------------------------------------------------------------------------
###Figure2H Code###
###TF Activity Analysis using Homer###
#First get differentially methylated CpGs from comparing S1 to S3 and S3 to S5
phenotype <- pheno(NK_RnB)
phenotype$Comparisons <- c(rep("all_else", 21),
                           rep("S1",4),
                           rep("all_else", 8),
                           rep("S3",4),
                           rep("all_else", 8),
                           rep("S5", 3)
)
slot(NK_RnB, "pheno") <- phenotype
num.cores <- 24 # adjust based on cores requested for HPC session
parallel.setup(num.cores)
dm <- rnb.execute.computeDiffMeth(NK_RnB, pheno.cols =  "Comparisons", 
                                  pheno.cols.all.pairwise = "Comparisons")
save.rnb.diffmeth(dm, "~/DNAmethylation/NKDI_NK_IV_Analysis/diff_meth")

#Find the probes that are hyper or hypomethylated between S1 vs S3 or S3 vs S5
dm <- load.rnb.diffmeth("~/DNAmethylation/NKDI_NK_IV_Analysis/diff_meth")
dmt_1 <- get.table(dm, get.comparisons(dm)[4], "sites", return.data.frame = T)
dmt_2 <- get.table(dm, get.comparisons(dm)[6], "sites", return.data.frame = T)
rownames(dmt_1) <- rownames(NK_RnB@sites)
rownames(dmt_2) <- rownames(NK_RnB@sites)
probes_1 <- na.omit(dmt_1[abs(dmt_1$mean.diff) >= 0.3, ])
probes_2 <- na.omit(dmt_2[abs(dmt_2$mean.diff) >= 0.3, ])
combined_probes <- unique(c(rownames(probes_1), rownames(probes_2)))
pref_probes <- read.csv(
  "~/DNAmethylation/NKDI_NK_IV_Analysis/805KAnnotation.csv")
combined_probes <- combined_probes[combined_probes %in% pref_probes$IlmnID]

#Find probe annotation data for 16559 differentially methylated probes
probe_anno <- annotation(NK_RnB, type = "sites", add.names = T)
probes_dm <- probe_anno[rownames(probe_anno) %in% combined_probes, ]

# Create BED file with ±100 bp windows
probes_bed <- data.frame(
  name = rownames(probes_dm),
  chrom = probes_dm$Chromosome,
  start = probes_dm$Start - 100,
  end = probes_dm$Start + 100,
  strand = probes_dm$Strand)
write.table(probes_bed, 
            "~/DNAmethylation/NKDI_NK_IV_Analysis/TF_Activity_Homer/16k_probes_100.txt", 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
#PERFORM HOMER ANALYSIS USING SLURM JOB ON HPC

# Read HOMER output correctly (using data.table for large files)
homer <- fread(
  "~/DNAmethylation/NKDI_NK_IV_Analysis/TF_Activity_Homer/New_Analysis/FK3vFK315/16k_motifs.txt", 
  skip = "PeakID (cmd=annotatePeaks.pl 16k_probes_100.txt hg38 -m ./known.motifs -cpu 40)", 
  header = TRUE)
colnames(homer)[1] <- "PeakID"

# Identify motif columns (those containing "/Homer Distance From Peak")
motif_cols <- grep("/Homer Distance From Peak", colnames(homer), value = TRUE)

# Extract clean motif names (everything before the first "/")
motif_names <- sub("^([^/]+).*", "\\1", motif_cols)

# Create mapping table of probes to motifs
motif_map <- rbindlist(
  lapply(seq_along(motif_cols), function(i) {
    col <- motif_cols[i]
    motif <- motif_names[i]
    # Get probes with this motif (non-empty distance values)
    probes <- homer[!is.na(get(col)) & get(col) != "", .(PeakID)]
    if(nrow(probes) > 0) data.table(PeakID = probes$PeakID, Motif = motif)
  }
  )
)
# Filter motifs with sufficient coverage (min 5 probes)
motif_counts <- motif_map[, .N, by = Motif]
keep_motifs <- motif_counts[N >= 5, Motif]
motif_map <- motif_map[Motif %in% keep_motifs]

# Verify
print(paste("Found", length(unique(motif_map$Motif)), "motifs with ≥5 probes"))
head(motif_map, 3)

# Load methylation data
beta_matrix <- meth(NK_RnB, row.names = TRUE)
colnames(beta_matrix) <- pheno(NK_RnB)$Sample_ID

#Find the probes that are hyper or hypomethylated between S1 vs S3 or S3 vs S5
dm <- load.rnb.diffmeth("~/DNAmethylation/NKDI_NK_IV_Analysis/diff_meth")
dmt_1 <- get.table(dm, get.comparisons(dm)[4], "sites", return.data.frame = T)
dmt_2 <- get.table(dm, get.comparisons(dm)[6], "sites", return.data.frame = T)
rownames(dmt_1) <- rownames(NK_RnB@sites)
rownames(dmt_2) <- rownames(NK_RnB@sites)
probes_1 <- na.omit(dmt_1[abs(dmt_1$mean.diff) >= 0.3, ])
probes_2 <- na.omit(dmt_2[abs(dmt_2$mean.diff) >= 0.3, ])
combined_probes <- unique(c(rownames(probes_1), rownames(probes_2)))
pref_probes <- read.csv("~/DNAmethylation/NKDI_NK_IV_Analysis/805KAnnotation.csv")
combined_probes <- combined_probes[combined_probes %in% pref_probes$IlmnID]

beta_dm <- beta_matrix[combined_probes, ]

# Calculate motif means per sample
# Get unique motifs from the mapping
unique_motifs <- unique(motif_map$Motif)

# Initialize sample x motif matrix
sample_motif_means <- matrix(NA,
                             nrow = ncol(beta_dm),
                             ncol = length(unique_motifs),
                             dimnames = list(colnames(beta_dm), unique_motifs))

# Fill matrix with mean methylation per motif per sample
for(m in unique_motifs) {
  probes <- motif_map[Motif == m, PeakID]
  valid_probes <- intersect(probes, rownames(beta_dm))
  if(length(valid_probes) > 0) {
    sample_motif_means[, m] <- colMeans(beta_dm[valid_probes, , drop = FALSE], na.rm = TRUE)
  }
}

# Calculate global methylation per sample (mean of all 16k probes)
universe_means <- colMeans(beta_dm, na.rm = TRUE)

# Probe-set normalize motif values
sample_motif_norm <- sample_motif_means / universe_means
rownames(sample_motif_norm) <- phenotype$Sample_ID

# If primary cells are available (e.g., for batch correction)
primary_celltypes <- c("S1","S2A","S2B","S3","S4A","S4B","S5")
primary_samples <- phenotype %>%
  filter(CellType %in% primary_celltypes) %>%
  pull(Sample_ID) %>%
  intersect(rownames(sample_motif_norm))

# Column-wise means (average per MOTIF across primary samples)
primary_avg <- colMeans(sample_motif_norm[primary_samples, ], na.rm = TRUE)

sample_tf_activity <- sweep(sample_motif_norm, 2, primary_avg, "-")

# Invert all values in the matrix/data.frame
sample_tf_activity_inverted <- -1 * sample_tf_activity

# Create metadata dataframe from rownames
metadata <- data.frame(
  Sample_ID = rownames(sample_tf_activity),
  CellType = phenotype$CellType
)

# Verify
metadata$CellType <- phenotype$CellType
head(metadata)

# Convert to long format and merge with metadata
tf_activity_long <- sample_tf_activity_inverted %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("Sample_ID") %>% 
  pivot_longer(
    cols = -Sample_ID,
    names_to = "Motif",
    values_to = "TF_Activity"
  ) %>% 
  left_join(metadata, by = "Sample_ID")

#Test if each motif’s activity in a cell type differs significantly from 0 (null = no activity change):
stat_results <- tf_activity_long %>% 
  group_by(CellType, Motif) %>% 
  summarise(
    Mean_Activity = mean(TF_Activity, na.rm = TRUE),
    p_value = t.test(TF_Activity, mu = 0)$p.value,
    .groups = "drop"
  ) %>% 
  mutate(p_adj = p.adjust(p_value, method = "fdr"))  # Correct for multiple testing

#Save a csv for importing the data into graphpad for normalizing and plotting
write.csv(tf_activity_long, 
          "~/DNAmethylation/NKDI_NK_IV_Analysis/tf_activity_D28.csv")
#------------------------------------------------------------------------------
###Figure 2I, Figure 5G and Figure 5H code###
###KIR Promoter Methylation Status###
#Look at KIR loci methylation at day 14 and 28
betas <- meth(NK_RnB, row.names=T)
colnames(betas) <- phenotype$CellType
betas <- as.data.frame(betas)

probes <- annotation(NK_RnB, type = "sites", add.names = T, include.regions = T)
genes <- annotation(NK_RnB, type = "genes", add.names = T)
promoters <- annotation(NK_RnB, type = 'promoters', add.names = T)

kir_gene_rows <- which(grepl("^KIR\\d", genes$symbol))
kir_probes_df <- probes[probes$genes %in% kir_gene_rows, ]
kir_promoter_rows <- which(grepl("^KIR\\d", promoters$symbol))
kir_promoter_probes_df <- probes[probes$promoters %in% kir_promoter_rows, ]

kir_probe_ids <- rownames(kir_probes_df)
kir_gene_betas <- betas[rownames(betas) %in% kir_probe_ids, ]
kir_promoter_ids <- rownames(kir_promoter_probes_df)
kir_promoter_betas <- betas[rownames(betas) %in% kir_promoter_ids, ]

avg_kir_methylation <- colMeans(kir_gene_betas)
avg_kir_promoter_methylation <- colMeans(kir_promoter_betas)

promoter_idx <- kir_promoter_probes_df[
  match(rownames(kir_promoter_betas),
        rownames(kir_promoter_probes_df)),
  "promoters"
]
gene_names <- ifelse(
  promoter_idx > 0,
  promoters$symbol[promoter_idx],
  NA
)
kir_promoter_betas_out <- cbind(
  Probe = rownames(kir_promoter_betas),
  Gene = gene_names,
  kir_promoter_betas
)

#Write a csv file for 
write.csv(
  kir_promoter_betas_out,
  "~/DNAmethylation/NKDI_NK_IV_Analysis/KIR_promoter_beta_values.csv",
  row.names = FALSE
)
#-----------------------------------------------------------------------------
###Figure 5F Code###
###Homer De Novo Motif Prediction at Day 14 between 15C and 15D###














