library(reshape2)
library(ggplot2)
library(tidyr)
library(dplyr)
library(broom)
############################
load("/Users/jeremyhappe/Desktop/Condorcet/TFE/syg.RData")
tree <- read.tree(file = "/Users/jeremyhappe/Desktop/Condorcet/TFE/BAMM_new_tree/full_tree.tre")
#############
Y <- syg
#Depending on what you want to select
Y <- na.omit(Y[, c("ID", "Species","Egg_size")])

common_taxa <- intersect(tree$tip.label, Y$Species)
if (length(common_taxa) == 0) {
  stop("No common taxa found between the tree and Excel file.")}
tree_pruned <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% common_taxa])
row.names(Y) <- Y$Species
Y <- Y[order(match(Y$Species, tree_pruned$tip.label)),]

#######################################
ln_DR <- (setNames((Y$Egg_size), rownames(Y)))
DR_contMap <- contMap(tree_pruned, ln_DR, res =1000, lwd = 2, type = "phylogram",outline=FALSE,fsize = 0.4)
Y$BAMM <- Y$speciation_rate #rename speciation rate
Y_subset <- Y[, c("DR", "BAMM", "ClaDs", "MISSE")]


cor_matrix <- cor(Y_subset, use = "complete.obs", method = "spearman") # spearman since it's not normaly distributed
cor_melt <- melt(cor_matrix) # Melt the correlation matrix into a long format

# Create the heatmap using ggplot2
ggplot(data = cor_melt, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "black") +  # Add correlation values
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.5, 
                       limit = c(0, 1), space = "Lab", name = "Correlation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(title = "Heatmap of Spearman Correlation Matrix", x = "", y = "")


Y_long <- Y %>%
  pivot_longer(cols = c(DR, ClaDs, MISSE), names_to = "Method", values_to = "Value")

# Calculate R² values for each method
models <- Y_long %>%
  group_by(Method) %>%
  do(model = lm(Value ~ BAMM, data = .))

# Extract R² values
model_summaries <- models %>%
  summarise(Method, r.squared = summary(model)$r.squared)

# Merge R² values with long-format data
Y_long <- left_join(Y_long, model_summaries, by = "Method")

# Create the faceted scatter plots
ggplot(Y_long, aes(x = BAMM, y = Value)) +
  geom_point(aes(color = Family)) +  # Adds the points
  geom_smooth(method = "lm", col = "blue", fill = "lightblue", level = 0.95) +  # Adds a linear regression line with confidence interval
  facet_wrap(~ Method, scales = "free_y") +  # Facet by the Method column
  labs(title = "Faceted Scatter Plots of DR vs Other Methods",
       x = "BAMM",
       y = "Value") +
  theme_minimal() +
  geom_text(data = model_summaries, aes(x = Inf, y = Inf, label = paste("R² =", round(r.squared, 2))),
            hjust = 1.1, vjust = 1.1, inherit.aes = FALSE)  # Add R² values to the plots


#### PLOT DR
library(phytools)

syg$Species <- syg$Species[match(tree$tip.label, syg$Species)]
ln_DR <- (setNames(log(syg$ClaDs), syg$Species))
DR_contMap <- contMap(tree, ln_DR, res =800, lwd = 3, type = "phylogram",outline=FALSE,)

setwd("~/Desktop/Condorcet/TFE/PLOT")
# Open a PNG device with increased size
png("DR_tree.png")  # res is the resolution in DPI
plot(DR_contMap,fsize=c(0.7,0.8),leg.txt="Diversification rate", lwd = 1, type = "phylogram", outline=FALSE)
dev.off()
