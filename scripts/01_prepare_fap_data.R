library(Seurat)
library(dplyr)
library(stringr)

# Load and subset the FAP object
fap <- readRDS("/Users/vayun/Downloads/SMG_delt_FAP.rds")
Idents(fap) <- "adjusted_celltype"
selected_cell_types <- c("CD74", "DLK1", "ATF3", "CD55", "GLI1")
fap_subset <- subset(fap, idents = selected_cell_types)

# Extract age from the sample ID
fap_subset$age <- as.numeric(str_extract(fap_subset$orig.ident, "(?<=_)[0-9]+"))

# Assign cell-cycle phases
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
fap_subset <- CellCycleScoring(fap_subset, s.features = s.genes,
                               g2m.features = g2m.genes, set.ident = FALSE)

# Calculate proliferative fraction per patient/sample
meta <- fap_subset[[]] %>% select(orig.ident, Phase)
prolif_df <- meta %>%
  mutate(Prolif = Phase %in% c("S", "G2M")) %>%
  group_by(orig.ident) %>%
  summarise(Prolif_Fraction = mean(Prolif))

# Add proliferative fraction back to the object
fap_subset$Prolif_Fraction <- prolif_df$Prolif_Fraction[match(fap_subset$orig.ident,
                                                              prolif_df$orig.ident)]

# (Optional) biological age formula from the repository
fap_subset$BioAge <- 35 - 100 * fap_subset$Prolif_Fraction

# Save subsetted object for downstream analysis
saveRDS(fap_subset, file = "data/fap_subset.rds")
