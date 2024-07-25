###################
### PPA example ###
###################
library(graph)
library(phylopath)
library(tidyverse)
library(dbplyr)
library(ape)
library(readxl)
library(DRstat)
library(TreeSimGM)
library(epm)
library(phytools)
library(caper)
library(ape)
library("expm")
library(readxl)
library(dplyr)
library(MuMIn)
library(ggplot2)
library(ggpubr)
load("/Users/jeremyhappe/Desktop/Condorcet/TFE/syg.RData")
tree <- read.tree(file = "/Users/jeremyhappe/Desktop/Condorcet/TFE/BAMM_new_tree/full_tree.tre")
########################################################
Y <- syg

#Select traits colour
Y <- na.omit(syg[, c(#"Polymorphic",
                      "Species", "DR", "speciation_rate","ClaDs","MISSE",
                       "Conspicuous",
                     #  "Red", 
                      # "Camouflage",
                     #  "Mimicry",
                       "Dots",
                       "Strips",
                     #  "Length",
                     'Stealth',
                     "Pattern"
                     
)])

#### SEX traits
Y <- na.omit(syg[, c("ID","Family", "Species" , "Egg_size", "Monogamous","speciation_rate","MISSE","DR","ClaDs",
                       "Length")])

common_taxa <- intersect(tree$tip.label, Y$Species)
if (length(common_taxa) == 0) stop("No common taxa found between the tree and Excel file.")
tree_pruned <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% common_taxa])

###FACTOR

Y$Monogamous <- ifelse(Y$Monogamous == 1, "Monogamous", "NO_Monogamous")
Y$Conspicuous <- ifelse(Y$Conspicuous == 1, "Conspicuous", "NO_Conspicuous")
Y$Camouflage <- ifelse(Y$Camouflage == 1, "Camouflage", "NO_Camouflage")
Y$Polymorphic <- ifelse(Y$Polymorphic == 1, "Polymorphic", "NO_Polymorphic")
Y$Red <- ifelse(Y$Red == 1, "Red", "NO_Red")
Y$Dots <- ifelse(Y$Dots == 1, "Dots", "N0_Dots")
Y$Strips <- ifelse(Y$Strips == 1, "Strips", "NO_Strips")
Y$Mimicry <- ifelse(Y$Mimicry == 1, "Mimicry", "NO_Mimicry")
Y$Pattern <- ifelse(Y$Pattern == 1, "Pattern", "NO_Pattern")
Y$Stealth <- ifelse(Y$Stealth == 1, "Stealth", "NO_Stealth")
Y$Sexual_dimorphism<- ifelse(Y$Sexual_dimorphism == 1, "SExu", "NO sex")
Y$Monogamous <- as.factor(Y$Monogamous)
Y$Conspicuous <- as.factor(Y$Conspicuous)
Y$Camouflage <- as.factor(Y$Camouflage)
Y$Polymorphic <- as.factor(Y$Polymorphic)
Y$Red<- as.factor(Y$Red)
Y$Dots<- as.factor(Y$Dots)
Y$Strips<- as.factor(Y$Strips)
Y$Mimicry<- as.factor(Y$Mimicry)
Y$Stealth<- as.factor(Y$Stealth)
Y$Pattern<- as.factor(Y$Pattern)
Y$Sexual_dimorphism <- as.factor(Y$Sexual_dimorphism)
str(Y)

# Current column names mapping to new names
new_names <- c(ID = "ID",
               Species ="Species",
               Family = "Family",
               Polymorphic = "P",
               Conspicuous = "G", # Make sure to correct the spelling if it's meant to be "Conspicuous"
               Red = "R",
               Camouflage = "Ca",
               Mimicry = "Mi",
               Dots = "D",
               Strips = "S",
               Length = "L",
               DR = "DR",
               Clutch_size = "C",
               Monogamous = "M",
               Pattern = "Z",
               Stealth = "St",
               speciation_rate = "BAMM",
               ClaDs = "ClaDS",
               MISSE = "MiSSE",
               Egg_size ="E",
               Brood_size = "B",
               Sexual_dimorphism ="X")

colnames(Y) <- new_names[colnames(Y)]

models2 <- define_model_set(
  Direct = c(ClaDS ~ G+S),
  Indirect = c(ClaDS~ G+S+Z),
  Pattern = c(ClaDS~G+Z),
 Pat_ind = c(ClaDS ~ G+Z+S,G~Z+St),
 Pat_ind2 = c(ClaDS ~ G+S,G~Z+St),
 
 )

model_clutch <- define_model_set(
  All = c( DR~E+M+L, E~M+L),
  indirect = c( DR~E+M+L),
  mid = c(DR~E+L+M, E~L)

)

# run ppa
data_frame <- as.data.frame(Y)
rownames(data_frame) <- data_frame$Species

result <- phylo_path(models2,data_frame, tree_pruned)

result <- phylo_path(model_clutch, data_frame, tree_pruned)
s <- summary(result)
plot(s)
b <- best(result)
plot(b)

# look at standardized coefficients and errors of the paths to look at a fitted model
coef_plot(b, error_bar = "se", order_by = "strength", to = "DR") +
  ggplot2::coord_flip()+
  ggplot2::theme_bw()








