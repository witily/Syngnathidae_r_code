library(readxl)
library("epm")
library("dplyr")
library("ggtree")

syg <- read_excel("/Users/jeremyhappe/Desktop/Condorcet/TFE/syngnathids_nodup.xlsx")
syg$Clutch_size <-as.numeric(syg$Clutch_size)

# Function to add underscores between words
add_underscores <- function(species_name) {
  gsub(" ", "_", species_name)
}

# Apply the function to the "Species" column
syg$Species <- sapply(syg$Species, add_underscores)

#####################
#GET BAMM

tree <- read.tree(file = "/Users/jeremyhappe/Desktop/Condorcet/TFE/BAMM_new_tree/full_tree.tre")
speciation_df <- read.table("/Users/jeremyhappe/Desktop/Condorcet/TFE/BAMM_new_tree/bamm_all_speciation_rate_df.txt", header = TRUE, sep = "\t")

# TREE CUT
common_taxa <- intersect(tree$tip.label, syg$Species)
syg <- syg %>% filter(Species %in% common_taxa)
syg <- syg[!duplicated(syg$Species), ]
syg <- merge(syg,speciation_df,  by.x = "Species", by.y = "ID", all.x=TRUE)

#GET DRSTAT

dr <- DRstat(tree) # Ensure the tree used here is the pruned tree
dr_df <- data.frame(ID = names(dr), DR = dr, row.names = NULL)
syg <- merge(syg, dr_df, by.x = "Species", by.y = "ID", all.x = TRUE)

#GET ClaDs

load("~/ClaDS.RData")
#syg$ClaDs <- CladsOutput$lambdatip_map
ClaDs <- CladsOutput$lambdatip_map
Species <- CladsOutput$tree$tip.label
ClaDs_df<-data.frame(Species,ClaDs)
syg <- merge(syg, ClaDs_df, by = "Species", all.x = TRUE)


#GET MISSE 

load("/Users/jeremyhappe/Desktop/Condorcet/TFE/MISSE_speciation.RData")

MISSE <- speciation_rates$speciation
Species <- speciation_rates$taxon
MISSE_df<-data.frame(Species,MISSE)
syg <- merge(syg, MISSE_df, by = "Species", all.x = TRUE)



#NEW values 

syg$Stealth <- as.integer(syg$Mimicry == 1 | syg$Camouflage == 1)
syg$Pattern <- as.integer(syg$Strips == 1| syg$Dots == 1)

#LOG transformation and shapiro test value

syg$DR<-log(syg$DR) # 0,0114
syg$ClaDs <- log(syg$ClaDs) #1.111e-06
syg$Length <- log(syg$Length) # 0,0006114
syg$Clutch_size <- log(syg$Clutch_size) # 0,00016
syg$Brood_size <- log(syg$Brood_size) ## 2.748e-05
syg$Egg_size <- log(syg$Egg_size) # 0.0255

# Facto transformation

cols_to_factor <- c("Family", "Monogamous", "Conspicuous", "Camouflage", "Polymorphic", 
                    "Red", "Dots", "Strips", "Mimicry", "Stealth", "Pattern", 
                    "Sexual_dimorphism")
syg[cols_to_factor] <- lapply(syg[cols_to_factor], as.factor)
str(syg)


save(syg, file = "/Users/jeremyhappe/Desktop/Condorcet/TFE/syg.RData")
load("/Users/jeremyhappe/Desktop/Condorcet/TFE/syg.RData")
########################################################

# Match the order of species in 'syg' to 'tree$tip.label'
syg <- syg[match(tree$tip.label, syg$Species), ]
# Create a named vector for renaming
rename_map <- c(
  "SY387" = "Eurypegasus_draconis",
  "SY025" = "Acentronura_cf._tentaculata",
  "SY328" = "Hippocampus_borboniensis"
)
# Rename species in 'tree$species' based on the map
syg_1$Species <- sapply(syg_1$ID, function(x) if (x %in% names(rename_map)) rename_map[x] else syg$Species[syg$ID == x])
# Verify the renaming
print(syg_1$Species)

# Set the tree's tip labels using the reordered 'Species' names
tree$tip.label <- syg_1$Species
# Verify the updated tree tip labels
print(tree$tip.label)









# Refresh

Y<- syg
Y<- as.data.frame(Y)
row.names(Y) <- Y$Species


#Y<- Y %>% filter(Family == "Syngnathidae")
Y <- na.omit(Y[, c("Species","Family","Length","DR","speciation_rate",
                   "Mimicry","Camouflage","Polymorphic","Red","Strips","Dots","Conspicuous","Appendages","Stealth","Pattern"
                   #"Monogamous",
                   #"Clutch_size" 
                   )])

common_taxa <- intersect(tree$tip.label, Y$Species)
if (length(common_taxa) == 0) {
  stop("No common taxa found between the tree and Excel file.")}
tree_pruned <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% common_taxa])
tree_pruned$node.label <- NULL



plot(tree)
#######################
###### PGLS
#######################
compdata <- comparative.data(tree_pruned, Y, 'Species', na.omit = FALSE)
bib<- pgls(DR ~ Clutch_size+ Monogamous, data = compdata, lambda = "ML")
summary(bib)
AIC(bib)

model_formulas <- list(
  DR ~ Polymorphic,
  DR ~ Polymorphic + Red,
  DR ~ Polymorphic + Length,
  DR ~ Polymorphic + Conspicuous + Strips,
  DR ~ Polymorphic + Conspicuous + Strips + Stealth,
  DR ~ Polymorphic + Length + Conspicuous,
  DR ~ Polymorphic + Length + Strips,
  DR ~ Polymorphic + Pattern + Stealth,
  DR ~ Polymorphic + Mimicry + Dots + Red
)

# Function to fit models and extract required information
extract_model_info <- function(formula, data) {
  model <- pgls(formula, data = data, lambda = "ML")
  summary_stats <- summary(model)
  aic_value <- AIC(model)
  # Assuming we are interested in the p-value of the first predictor
  first_pred_p_value <- ifelse(ncol(summary_stats$coefficients) > 3, summary_stats$coefficients[2, 4], NA)
  return(matrix(c(summary_stats$r.squared, summary_stats$adj.r.squared, aic_value, first_pred_p_value), nrow = 1))
}

# Apply the function to each model and create a data frame
model_info <- do.call(rbind, lapply(model_formulas, extract_model_info, data = compdata))
colnames(model_info) <- c("R-squared", "Adjusted R-squared", "AIC", "P-value (1st Predictor)")
rownames(model_info) <- c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6", "Model 7", "Model 8", "Model 9")

# Print the resulting data frame
print(model_info)

model_info <- as.data.frame(model_info)
model_info$Model <- row.names(model_info)
library(writexl)
write_xlsx(model_info, "/Users/jeremyhappe/Desktop/Condorcet/TFE/Model_AIC.xlsx")




