library(TreeSimGM)
library(epm)
library(phytools)
library(caper)
library(ape)
library("expm")
library(readxl)
library(dplyr)
library("phangorn")
library("diversitree")


load("/Users/jeremyhappe/Desktop/Condorcet/TFE/syg.RData")
tree <- read.tree(file = "/Users/jeremyhappe/Desktop/Condorcet/TFE/BAMM_new_tree/full_tree.tre")
########


Y <- na.omit(syg[, c("ID", "Species","DR","Family",
                     "Dots"
                     )])

common_taxa <- intersect(tree$tip.label, Y$Species)
if (length(common_taxa) == 0) {
  stop("No common taxa found between the tree and Excel file.")}
tree_pruned <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% common_taxa])
tree_pruned$node.label <- NULL
#######################################
#####END SET UP +++++ DELTA
############################################

tree_pruned$edge.length[tree_pruned$edge.length==0] <- quantile(tree_pruned$edge.length,0.1)*0.1

Y1 <- Y[order(match(Y$Species, tree_pruned$tip.label)),]
row.names(Y1) <- Y1$Species

fmode<-setNames(Y1$Dots,rownames(Y1))

###################################
#### Ancestral state
###################################


fmode<-setNames(Y1$Egg_size,rownames(Y1))

##### COMPARAISON OF"ER" 
######

anc_states <- fastAnc(tree_pruned, trait_data, vars=TRUE, CI=TRUE)
# Extract ancestral states and confidence intervals
reconstructed_states <- anc_states$ace
confidence_intervals <- anc_states$CI
print(reconstructed_states, printlen=10)
# Plot the phylogenetic tree
plot(tree_pruned, show.node.label=TRUE)



consp_map2<- make.simmap(tree_pruned, fmode,"ER", nsim=100)

map_summary <- describe.simmap(consp_map2)

## plot a posterior probabilities of ancestral states
cols<-setNames(c("blue","red"),levels(fmode))
# Plot the average of the summarized map with smaller pie charts and tip labels
plot(map_summary, colors = cols, ftype = "reg", fsize = 0.5, cex.pie = 0.9)

# Adjust plot margins to make space for the legend
par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)

legend("topright", inset = c(2, 0), # Use inset to move the legend inside the margin
       legend = c("Non-Stealth", "Stealth"), 
       pch = 21, pt.bg = cols, pt.cex = 1.5, cex = 0.8)

## plot posterior density on the number of changes
plot(density(consp_map2),bty="l")
title(main="Posterior distribution of changes of each type",
      font.main=3)
## End(Not run)



#### D stat

lambda<-phylosig(tree_pruned,Y1$Egg_size,method = "lambda",nsim = 100, test = TRUE)
quartz()

#### PAGELS 
a<-phylo.d(Y1, tree_pruned, Species, Dots, permut = 1000, rnd.bias=NULL)
comp_data <- comparative.data(tree_pruned, Y1, Species, vcv = TRUE)
lambda_result <- pgls(Egg_size ~ 1, comp_data, lambda = 'ML')
summary(lambda_result)


