library(caper)
library(dplyr)
library(readxl)
library(ape)
####################
load("/Users/jeremyhappe/Desktop/Condorcet/TFE/syg.RData")
tree <- read.tree(file = "/Users/jeremyhappe/Desktop/Condorcet/TFE/BAMM_new_tree/full_tree.tre")

###################
###################

# Define binary traits to be analyzed

traits <- c("Monogamous", "Red", "Conspicuous", "Strips",  
            "Mimicry", "Camouflage", "Stealth", "Polymorphic", "Pattern", 
            "Sexual_dimorphism")


# Function to perform PGLS for a given trait
perform_pgls <- function(trait, data, tree) {
  # Select the relevant columns
  selected_data <- na.omit(data[, c("Species", "ClaDs", trait)])
  
  # Prune the tree for the selected species
  common_taxa <- intersect(tree$tip.label, selected_data$Species)
  if (length(common_taxa) == 0) {
    stop(paste("No common taxa found between the tree and the data for trait:", trait))
  }
  
  tree_pruned <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% common_taxa])
  tree_pruned$node.label <- NULL
  
  # Create comparative data object
  compdata <- comparative.data(tree_pruned, selected_data, 'Species', na.omit = FALSE)
  
  # Perform PGLS analysis
  formula <- as.formula(paste("ClaDs ~", trait))
  model <- pgls(formula, data = compdata, lambda ="ML")
  
  # Print the summary of the model
  cat("\nSummary for trait:", trait, "\n")
  print(summary(model))
}

# Loop through each trait and perform PGLS
for (trait in traits) {
  perform_pgls(trait, syg, tree)
}

#####################
#####################

# Define continuous traits to be analyzed
traits <- c("Clutch_size","Length", "Egg_size", "Brood_size")

# Function to perform PGLS for a given trait
perform_pgls <- function(trait, data, tree) {
  # Select the relevant columns
  selected_data <- na.omit(data[, c("Species", "ClaDs", trait)])
  
  # Prune the tree for the selected species
  common_taxa <- intersect(tree$tip.label, selected_data$Species)
  if (length(common_taxa) == 0) {
    stop(paste("No common taxa found between the tree and the data for trait:", trait))
  }
  
  tree_pruned <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% common_taxa])
  tree_pruned$node.label <- NULL
  
  # Create comparative data object
  compdata <- comparative.data(tree_pruned, selected_data, 'Species', na.omit = FALSE)
  
  # Perform PGLS analysis
  formula <- as.formula(paste("ClaDs ~", trait))
  model <- pgls(formula, data = compdata, lambda = "ML")
  
  # Print the summary of the model
  cat("\nSummary for trait:", trait, "\n")
  print(summary(model))
}

# Loop through each trait and perform PGLS
for (trait in traits) {
  perform_pgls(trait, syg, tree)
}


######################
#### CUMULATIVE models
######################

Y <- na.omit(syg[, c("Species","Family","Length","DR","speciation_rate","ClaDs","MISSE",
                    "Mimicry","Camouflage","Red","Strips","Dots","Conspicuous","Stealth","Pattern"
                    ,"Polymorphic"
                #   "Monogamous","Egg_size","Length"
                
)])
   Y<-na.omit(syg[, c("Species","Family","Dots","DR","speciation_rate","ClaDs","MISSE"
                  
                      #   "Monogamous","Egg_size","Length"
                      
   )])
#Prune_tree
common_taxa <- intersect(tree$tip.label, Y$Species)
if (length(common_taxa) == 0) {
  stop("No common taxa found between the tree and Excel file.")}
tree_pruned <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% common_taxa])
tree_pruned$node.label <- NULL

#Create a compdata
compdata <- comparative.data(tree_pruned, Y, 'Species', na.omit = FALSE)
compdata$data$Pattern <- scale(compdata$data$Pattern)
compdata$data$ClaDs <- scale(compdata$data$ClaDs)
compdata$data$ClaDs <- (compdata$data$ClaDs)/2
compdata$data$Length <- (compdata$data$Length)/2
compdata$data$Length <- scale(compdata$data$Length)
a<- pgls(speciation_rate ~ Dots, data = compdata, lambda = "ML")
summary (a)
#define the cumulative models based on the PGLS results
a
a <- pgls(ClaDs ~ Strips+ Length+ Camouflage+ Mimicry+ Polymorphic + Red + Dots+ Conspicuous+ Pattern, data=compdata, lambda = "ML")
a<- pgls(ClaDs ~ Pattern +  Stealth + Length + Conspicuous, data = compdata,, lambda = 'ML')
a<- pgls(ClaDs ~ 1, data = compdata)
summary (a)
a<- pgls(ClaDs~1, data = compdata, lambda = "ML")
summary (a)
a<- pgls(ClaDs ~ Length, data = compdata, lambda = "ML")
summary (a)

a$phyres
AIC(a)

#ADVICE : create models seperatly for a better global view 

model_formulas <- list(
  DR ~ Pattern,
  DR ~ Pattern + Conspicuous + Mimicry +Polymorphic ,
  DR ~ Dots+ Conspicuous + Mimicry,
  DR ~ Strips +Conspicuous + Mimicry,
  DR ~  Pattern+ Conspicuous + Camouflage ,
  DR ~  Pattern +  Conspicuous + Mimicry,
  DR ~  Pattern +  Stealth + Conspicuous ,
  DR ~  Red,
  DR ~ Strips 
)

model_formulas <- list(
  speciation_rate ~ Mimicry,
  speciation_rate ~ Pattern+ Conspicuous,
  speciation_rate ~ Conspicuous ,
  speciation_rate ~ Camouflage,
  speciation_rate ~ Length,
  speciation_rate ~ Strips ,
  speciation_rate ~ Length + Strips,
  speciation_rate ~ Stealth ,
  speciation_rate ~ Dots+ Strips
)

model_formulas <- list(
  ClaDs ~ Mimicry,
  ClaDs ~ Pattern,
  ClaDs ~ Conspicuous ,
  ClaDs ~ Camouflage,
  ClaDs ~ Strips ,
  ClaDs ~ Red,
  ClaDs ~ Stealth ,
  ClaDs ~ Dots
)
a<- pgls(speciation ~ Stealth, data = compdata)

summary (a)
AIC(a)

a<- pgls(DR ~ Pattern +  Stealth + Conspicuous + Mimicry
, data = compdata, lambda = "ML" )
a<- pgls(DR ~ Pattern + Conspicuous , data = compdata, lambda = "ML" )
a<- pgls(DR ~ Pattern + Conspicuous + Mimicry , data = compdata, lambda = "ML" )
a<- pgls(DR ~ Pattern + Conspicuous + Camouflage, data = compdata, lambda = "ML" )
a<- pgls(speciation_rate ~ Length+ Camouflage+ Mimicry+ Polymorphic + Red + Pattern+ Conspicuous, data = compdata, lambda = "ML" )

# Function to fit models and extract required information

############# loop for binary variables 

calculate_AIC <- function(model) {
  return(AIC(model))
}

# Initial empty model //// change DR accordingly with BAMM /CLADS/MISSE
current_model <- pgls(DR ~ 1, data = compdata)
best_AIC <- calculate_AIC(current_model)

# List of predictors
predictors <- c("Mimicry", "Camouflage", "Red", "Strips", "Dots", "Conspicuous", "Stealth", "Pattern")
selected_predictors <- c()

# List to store all candidate models and their AICs
all_models <- list(current_model)
all_AICs <- c(best_AIC)

# Forward stepwise selection
for (i in 1:length(predictors)) {
  candidate_models <- list()
  candidate_AICs <- c()
  
  for (predictor in predictors) {
    if (!(predictor %in% selected_predictors)) {
      formula <- as.formula(paste("DR ~", paste(c(selected_predictors, predictor), collapse = "+")))
      candidate_model <- pgls(formula, data = compdata)
      candidate_models[[predictor]] <- candidate_model
      candidate_AICs <- c(candidate_AICs, calculate_AIC(candidate_model))
    }
  }
  
  # Output the AICs of the candidate models at each step
  cat("\nStep", i, "candidate models' AICs:\n")
  print(candidate_AICs)
  
  # Select the best model from candidates
  best_candidate_index <- which.min(candidate_AICs)
  if (candidate_AICs[best_candidate_index] < best_AIC) {
    best_AIC <- candidate_AICs[best_candidate_index]
    best_predictor <- names(candidate_models)[best_candidate_index]
    selected_predictors <- c(selected_predictors, best_predictor)
    current_model <- candidate_models[[best_predictor]]
  }
  
  # Store all candidate models and their AICs
  all_models <- c(all_models, candidate_models)
  all_AICs <- c(all_AICs, candidate_AICs)
}

# Sort the models by AIC and keep the 3 best models
sorted_indices <- order(all_AICs)
all_models <- all_models[sorted_indices]
all_AICs <- all_AICs[sorted_indices]

# Remove duplicate models
unique_formulas <- unique(sapply(all_models, function(model) deparse(formula(model))))
unique_indices <- match(unique_formulas, sapply(all_models, function(model) deparse(formula(model))))
all_models <- all_models[unique_indices]
all_AICs <- all_AICs[unique_indices]

# Select the top 3 unique models
top_3_models <- all_models[1:min(3, length(all_models))]
top_3_AICs <- all_AICs[1:min(3, length(all_AICs))]

# Display summaries of the top 3 models
for (i in 1:length(top_3_models)) {
  cat("\nModel", i, "with AIC:", top_3_AICs[i], "\n")
  print(summary(top_3_models[[i]]))
}




#################
#################

############# loop for binary variables 

calculate_AIC <- function(model) {
  return(AIC(model))
}

# Initial empty model //// change DR accordingly with BAMM /CLADS/MISSE
current_model <- pgls(DR ~ 1, data = compdata)
best_AIC <- calculate_AIC(current_model)

# List of predictors
predictors <- c("")
selected_predictors <- c()

# List to store all candidate models and their AICs
all_models <- list(current_model)
all_AICs <- c(best_AIC)

# Forward stepwise selection
for (i in 1:length(predictors)) {
  candidate_models <- list()
  candidate_AICs <- c()
  
  for (predictor in predictors) {
    if (!(predictor %in% selected_predictors)) {
      formula <- as.formula(paste("DR ~", paste(c(selected_predictors, predictor), collapse = "+")))
      candidate_model <- pgls(formula, data = compdata)
      candidate_models[[predictor]] <- candidate_model
      candidate_AICs <- c(candidate_AICs, calculate_AIC(candidate_model))
    }
  }
  
  # Output the AICs of the candidate models at each step
  cat("\nStep", i, "candidate models' AICs:\n")
  print(candidate_AICs)
  
  # Select the best model from candidates
  best_candidate_index <- which.min(candidate_AICs)
  if (candidate_AICs[best_candidate_index] < best_AIC) {
    best_AIC <- candidate_AICs[best_candidate_index]
    best_predictor <- names(candidate_models)[best_candidate_index]
    selected_predictors <- c(selected_predictors, best_predictor)
    current_model <- candidate_models[[best_predictor]]
  }
  
  # Store all candidate models and their AICs
  all_models <- c(all_models, candidate_models)
  all_AICs <- c(all_AICs, candidate_AICs)
}

# Sort the models by AIC and keep the 3 best models
sorted_indices <- order(all_AICs)
all_models <- all_models[sorted_indices]
all_AICs <- all_AICs[sorted_indices]

# Remove duplicate models
unique_formulas <- unique(sapply(all_models, function(model) deparse(formula(model))))
unique_indices <- match(unique_formulas, sapply(all_models, function(model) deparse(formula(model))))
all_models <- all_models[unique_indices]
all_AICs <- all_AICs[unique_indices]

# Select the top 3 unique models
top_3_models <- all_models[1:min(3, length(all_models))]
top_3_AICs <- all_AICs[1:min(3, length(all_AICs))]

# Display summaries of the top 3 models
for (i in 1:length(top_3_models)) {
  cat("\nModel", i, "with AIC:", top_3_AICs[i], "\n")
  print(summary(top_3_models[[i]]))
}



a <- pgls(ClaDs ~ Brood_size, data = compdata ,lambda = 0.5043)
summary(a)
AIC(a)
