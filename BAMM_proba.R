
tree <- read.tree(file = "/Users/jeremyhappe/Desktop/Condorcet/TFE/full_tree.tre")
syg <- read_excel("/Users/jeremyhappe/Desktop/Condorcet/TFE/syngnathids_nodup.xlsx")

common_taxa <- intersect(tree$tip.label, syg$Species)
if (length(common_taxa) == 0) stop("No common taxa found between the tree and Excel file.")
tree_pruned <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% common_taxa])


Y <- syg %>% filter(ID %in% common_taxa)


Y <- Y[!duplicated(Y$Species), ]

common_taxa <- intersect(tree_pruned$tip.label, Y$ID)
if (length(common_taxa) == 0) {
  stop("No common taxa found between the tree and Excel file.")}
tree_pruned <- drop.tip(tree_pruned, tree_pruned$tip.label[!tree_pruned$tip.label %in% common_taxa])


### with the addition of some "theoretical new species" such as sp_2,.... i might have gone a bit too quick on that 
Y_use <- Y %>%
  mutate(samplingFraction = case_when(
    Family == "Syngnathidae" ~ (159/329),
    Family == "Solenostomidae" ~ (4/6),
    Family == "Pegasidae" ~ (2/8),
    Family == "Mullidae" ~ (51/107),
    Family == "Fistulariidae" ~ (4/5),
    Family == "Dactylopteridae" ~ (6/7),
    Family == "Callionymidae" ~ (33/201),
    Family == "Centriscidae" ~ (7/13),
    Family == "Aulostomidae" ~ (3/3),
    TRUE ~ NA_real_  # This line is for families not listed above, assigning NA
  ))


save_pruned_tree <- function(tree, filepath) {
  write.tree(tree, file = filepath)
}

Bamm_proba<- na.omit(Y_use[, c("Species","Family", "samplingFraction")])
row.names(Bamm_proba)<- Y_use$Species
Bamm_proba <- Bamm_proba %>% rename(cladeName = Family, speciesName = Species)                       

write.tree(tree_pruned,"/Users/jeremyhappe/Desktop/Condorcet/TFE/tree_pruned.tre" )                       
write_tsv(Bamm_proba,"/Users/jeremyhappe/Desktop/Condorcet/TFE/Bamm_proba.tsv")            
                       
                       
                       