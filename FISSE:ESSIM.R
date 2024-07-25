load("/Users/jeremyhappe/Desktop/Condorcet/TFE/syg.RData")
tree <- read.tree(file = "/Users/jeremyhappe/Desktop/Condorcet/TFE/BAMM_new_tree/full_tree.tre")
########
library(BAMMtools)
library(dplyr)
library(ape)
library(diversitree)

# Prepare trait data
Y <- na.omit(syg[, c("ID", "Species", "Egg_size")])

# Find common taxa and prune the tree
common_taxa <- intersect(tree$tip.label, Y$Species)
if (length(common_taxa) == 0) {
  stop("No common taxa found between the tree and Excel file.")
}
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, common_taxa))

# Ensure the tree is ultrametric
tree_ultrametric <- chronos(tree_pruned)

# Filter and arrange trait data to match the pruned tree
Y_filtered <- Y %>% filter(Species %in% common_taxa) %>% arrange(match(Species, tree_pruned$tip.label))

# Set row names for the trait data
row.names(Y_filtered) <- Y_filtered$Species

# Prepare binary trait data
fmode <- setNames(Y_filtered$Egg_size, row.names(Y_filtered))

# Run the FISSE.binary function
res <- FISSE.binary(tree_ultrametric, fmode)

row_name = "Sexual Dimorphism"

# Add new results to the dataframe with a row name
results_df <- add_results_to_df(res, results_df, row_name)

##################
##################

library("writexl")
results_df$Traits <- row.names(results_df)
write_xlsx(results_df, "/Users/jeremyhappe/Desktop/Condorcet/TFE/FISSE_TOT.xlsx")


##########
##.  ESSIM
#########
library("hisse")

# Prepare trait data
Y <- na.omit(syg[, c("ID", "Species", "Egg_size")])

# Find common taxa and prune the tree
common_taxa <- intersect(tree$tip.label, Y$Species)
if (length(common_taxa) == 0) {
  stop("No common taxa found between the tree and Excel file.")
}
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, common_taxa))
Y_filtered <- Y %>% filter(Species %in% common_taxa) %>% arrange(match(Species, tree_pruned$tip.label))
row.names(Y_filtered) <- Y_filtered$Species


fmode <- setNames(Y_filtered$Egg_size, row.names(Y_filtered))

# Run the essim function
essim(tree_pruned, fmode, nsim = 10000)
str(fmode)

a<-cor.test(syg$DR ,syg$speciation_rate)
a$p.value

