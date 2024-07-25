### HERE IS THE CODE TO MAKE THE CORRELATION HEATMAPS AND QQPLOT
###
library(corrplot)
library(Hmisc)

load("/Users/jeremyhappe/Desktop/Condorcet/TFE/syg.RData")
tree <- read.tree(file = "/Users/jeremyhappe/Desktop/Condorcet/TFE/BAMM_new_tree/full_tree.tre")

Y <- syg

scatter_plot <- ggplot(Y, aes(x = DR, y = ClaDs, color = Family)) +
  geom_point(size = 3) +
  labs(title = "DR vs ClaDS",
       x = "DR",
       y = "ClaDS") +
  theme_minimal()


# Compute the correlation between speciation_rate and ClaDS
cor_value <- cor(Y$ClaDs, Y$MISSE, use = "complete.obs")

# Create a correlation plot
cor_plot <- ggplot(Y, aes(x = ClaDs, y = MISSE)) +
  geom_point(aes(color = Family), size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  annotate("text", x = Inf, y = Inf, label = paste("Correlation:", round(cor_value, 2)),
           hjust = 1.1, vjust = 2, size = 5, color = "black") +
  labs(title = "Correlation Plot of Speciation Rate vs ClaDS",
       x = "ClaDS",
       y = "MISSE") +
  theme_minimal() +
  theme(legend.position = "right")

#Name transformation 
Y_subs <- syg[, c("DR", "speciation_rate", "ClaDs", "MISSE")]
Y_subset <- syg[, c("DR","MISSE")]

Y_subset$BAMM <- Y_subs$speciation_rate
Y_subset$ClaDS <- Y_subs$ClaDs


correlation_matrix <- rcorr(as.matrix(Y_subset), type="spearman")
R_values <- correlation_matrix$r
R2_values <- R_values^2

# Calculate R^2 values from linear models
R2_values <- matrix(NA, nrow=ncol(Y_subset), ncol=ncol(Y_subset))
rownames(R2_values) <- colnames(Y_subset)
colnames(R2_values) <- colnames(Y_subset)

for (i in 1:ncol(Y_subset)) {
  for (j in i:ncol(Y_subset)) {
    if (i != j) {
      model <- lm(Y_subset[, i] ~ Y_subset[, j])
      R2_values[i, j] <- summary(model)$r.squared
      R2_values[j, i] <- R2_values[i, j] # Make the matrix symmetric
    } else {
      R2_values[i, j] <- 1 # R^2 of a variable with itself is 1
    }
  }
}

# Combine correlation and R^2 into a single matrix for heatmap
combined_matrix <- matrix(paste0(
  "R: ", round(R_values, 2), "\n",
  "R^2: ", round(R2_values, 2)
), nrow=nrow(R_values), ncol=ncol(R_values))
rownames(combined_matrix) <- colnames(Y_subset)
colnames(combined_matrix) <- colnames(Y_subset)

col_palette <- colorRampPalette(c("white", "blue"))(200)

# Mask values below 0.5 for R_values
R_values[R_values < 0.5] <- NA

# Create a heatmap for R and R^2 values
corrplot(
  R_values, method="color", type="upper",
  addCoef.col = "black", # Add correlation coefficients
  tl.col="black", tl.srt=45, # Text label color and rotation
  title = "Spearman Correlation Coefficients and R^2 Values",
  col = col_palette, is.corr = FALSE, na.label = "NA"
)

# Extract p-values
p_values <- correlation_matrix$P

# Combine correlation and R^2 into a single matrix for heatmap
combined_matrix <- matrix(paste0(
  "R: ", round(R_values, 2), "\n",
  "R^2: ", round(R2_values, 2)
), nrow=nrow(R_values), ncol=ncol(R_values))
rownames(combined_matrix) <- colnames(Y_subset)
colnames(combined_matrix) <- colnames(Y_subset)

# Define color palette from 0.5 to 1
col_palette <- colorRampPalette(c("white", "Red"))(200)

# Create a heatmap for R and R^2 values
corrplot(
  R_values, method="color", type="upper",
  addCoef.col = "black", # Add correlation coefficients
  tl.col="black", tl.srt=45, # Text label color and rotation
  title = "Correlation Coefficients and R^2 Values",
  col = col_palette, is.corr = FALSE, na.label = "NA"
)

corrplot(
  R2_values, method="color", type="upper",
  addCoef.col = "black", # Add R^2 values
  tl.col="black", tl.srt=45, # Text label color and rotation
  title = "Heatmap of R^2 Values from Linear Models",
  col = col_palette, is.corr = FALSE, na.label = "NA"
)


# Fit linear models
model_bamm <- lm(DR ~ BAMM, data = Y_subset)
model_clads <- lm(DR ~ ClaDS, data = Y_subset)
model_misse <- lm(DR ~ MISSE, data = Y_subset)

# Plot diagnostic plots
par(mfrow = c(3, 2))
quartz()
plot(model_bamm, which = 1:2, main = "BAMM Model")
plot(model_clads, which = 1:2, main = "ClaDS Model")
plot(model_misse, which = 1:2, main = "MISSE Model")


model_bamm <- lm(DR ~ BAMM, data = Y_subset)
model_clads <- lm(DR ~ ClaDS, data = Y_subset)
model_misse <- lm(DR ~ MISSE, data = Y_subset)
model_bamm_clads <- lm(BAMM ~ ClaDS, data = Y_subset)
model_bamm_misse <- lm(BAMM ~ MISSE, data = Y_subset)
model_clads_misse <- lm(ClaDS ~ MISSE, data = Y_subset)



# BAMM Model
res_bamm <- resid(model_bamm)
fit_bamm <- fitted(model_bamm)

# ClaDS Model
res_clads <- resid(model_clads)
fit_clads <- fitted(model_clads)

# MISSE Model
res_misse <- resid(model_misse)
fit_misse <- fitted(model_misse)

# BAMM ~ ClaDS Model
res_bamm_clads <- resid(model_bamm_clads)
fit_bamm_clads <- fitted(model_bamm_clads)

# BAMM ~ MISSE Model
res_bamm_misse <- resid(model_bamm_misse)
fit_bamm_misse <- fitted(model_bamm_misse)

# ClaDS ~ MISSE Model
res_clads_misse <- resid(model_clads_misse)
fit_clads_misse <- fitted(model_clads_misse)

# Residuals vs Fitted for BAMM
p1 <- ggplot(data.frame(Fitted = fit_bamm, Residuals = res_bamm), aes(Fitted, Residuals)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  ggtitle("DR~BAMM: Residuals vs Fitted")

# Normal Q-Q for BAMM
p2 <- ggplot(data.frame(Sample = res_bamm), aes(sample = Sample)) +
  stat_qq() +
  stat_qq_line() +
  ggtitle("DR~BAMM: Normal Q-Q Plot")

# Residuals vs Fitted for ClaDS
p3 <- ggplot(data.frame(Fitted = fit_clads, Residuals = res_clads), aes(Fitted, Residuals)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  ggtitle("DR~ClaDS: Residuals vs Fitted")

# Normal Q-Q for ClaDS
p4 <- ggplot(data.frame(Sample = res_clads), aes(sample = Sample)) +
  stat_qq() +
  stat_qq_line() +
  ggtitle("DR~ClaDS: Normal Q-Q Plot")

# Residuals vs Fitted for MISSE
p5 <- ggplot(data.frame(Fitted = fit_misse, Residuals = res_misse), aes(Fitted, Residuals)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  ggtitle("DR~MISSE: Residuals vs Fitted")

# Normal Q-Q for MISSE
p6 <- ggplot(data.frame(Sample = res_misse), aes(sample = Sample)) +
  stat_qq() +
  stat_qq_line() +
  ggtitle("DR~MISSE: Normal Q-Q Plot")

# Residuals vs Fitted for BAMM ~ ClaDS
p7 <- ggplot(data.frame(Fitted = fit_bamm_clads, Residuals = res_bamm_clads), aes(Fitted, Residuals)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  ggtitle("BAMM~ClaDS: Residuals vs Fitted")

# Normal Q-Q for BAMM ~ ClaDS
p8 <- ggplot(data.frame(Sample = res_bamm_clads), aes(sample = Sample)) +
  stat_qq() +
  stat_qq_line() +
  ggtitle("BAMM~ClaDS: Normal Q-Q Plot")

# Residuals vs Fitted for BAMM ~ MISSE
p9 <- ggplot(data.frame(Fitted = fit_bamm_misse, Residuals = res_bamm_misse), aes(Fitted, Residuals)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  ggtitle("BAMM~MISSE: Residuals vs Fitted")

# Normal Q-Q for BAMM ~ MISSE
p10 <- ggplot(data.frame(Sample = res_bamm_misse), aes(sample = Sample)) +
  stat_qq() +
  stat_qq_line() +
  ggtitle("BAMM~MISSE: Normal Q-Q Plot")

# Residuals vs Fitted for ClaDS ~ MISSE
p11 <- ggplot(data.frame(Fitted = fit_clads_misse, Residuals = res_clads_misse), aes(Fitted, Residuals)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  ggtitle("ClaDS~MISSE: Residuals vs Fitted")

# Normal Q-Q for ClaDS ~ MISSE
p12 <- ggplot(data.frame(Sample = res_clads_misse), aes(sample = Sample)) +
  stat_qq() +
  stat_qq_line() +
  ggtitle("ClaDS~MISSE: Normal Q-Q Plot")
library(gridExtra)
quartz()
# Arrange the plots in a facet-like structure
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, ncol = 4)


########################################
par(mfrow = c(2, 2)) # Set up the plotting area to have 4 plots (2x2)

# QQ plot for DR
qqnorm(Y$DR, main = "QQ Plot of DR")
qqline(Y$DR, col = "red")

# QQ plot for BAMM
qqnorm(Y$speciation_rate, main = "QQ Plot of BAMM")
qqline(Y$speciation_rate, col = "red")

# QQ plot for ClaDs
qqnorm(Y$ClaDs, main = "QQ Plot of ClaDs")
qqline(Y$ClaDs, col = "red")

# QQ plot for MISSE
qqnorm(Y$MISSE, main = "QQ Plot of MISSE")
qqline(Y$MISSE, col = "red")

######################################## 

Y<- syg
Y$BAMM <- Y$speciation_rate

Y_long <- Y %>%
  pivot_longer(cols = c(DR,BAMM, ClaDs, MISSE), names_to = "Method", values_to = "Value")

# Create the facet plot
ggplot(Y_long, aes(y = Family, x = Value, fill = Family)) +
  geom_boxplot() +
  geom_point(alpha = 0.3, color = "black", width = 0.2) +
  facet_wrap(~ Method, scales = "free_x") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1)) +
  labs(title = "Comparison of BAMM, MISSE, and ClaDs across Families",
       y = "Family",
       x = "Value")

####################

#Prune_tree
common_taxa <- intersect(tree$tip.label, Y$Species)
if (length(common_taxa) == 0) {
  stop("No common taxa found between the tree and Excel file.")}
tree_pruned <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% common_taxa])
tree_pruned$node.label <- NULL

#Create a compdata
compdata <- comparative.data(tree_pruned, Y, 'Species', na.omit = FALSE)

a<-pgls(ClaDs~MISSE, data= compdata, lambda ="ML")
summary(a)
# Extract fitted values and residuals
fitted_values_a <- fitted(a)
residuals_a <- residuals(a)

# Create residuals vs. fitted values plot
plot(fitted_values_a, residuals_a,
     xlab = "Fitted Values",
     ylab = "Residuals",
     main = "Residuals vs Fitted Values",
     pch = 20, col = "blue")
abline(h = 0, col = "red", lwd = 2)
library(car) 
# Create QQ plot
quartz()
qqPlot(residuals_a, main = "QQ Plot of Residuals for PGLS Model", ylab = "Residuals")

Y_subs <- Y[, c("DR", "BAMM", "ClaDs", "MISSE")]
Y_subset <- Y[, c("DR", "MISSE")]

Y_subset$BAMM <- Y_subs$BAMM
Y_subset$ClaDs <- Y_subs$ClaDs

# Create correlation matrix
correlation_matrix <- rcorr(as.matrix(Y_subset), type="spearman")
R_values <- correlation_matrix$r
R2_values <- R_values^2

# Create a matrix to store p-values and adjusted R^2 values
pgls_p_values <- matrix(NA, nrow=ncol(Y_subset), ncol=ncol(Y_subset))
pgls_adj_R2_values <- matrix(NA, nrow=ncol(Y_subset), ncol=ncol(Y_subset))
rownames(pgls_p_values) <- colnames(Y_subset)
colnames(pgls_p_values) <- colnames(Y_subset)
rownames(pgls_adj_R2_values) <- colnames(Y_subset)
colnames(pgls_adj_R2_values) <- colnames(Y_subset)


# Loop through each pair of variables and perform PGLS
for (i in 1:ncol(Y_subset)) {
  for (j in 1:ncol(Y_subset)) {
    if (i != j) {
      formula <- as.formula(paste(colnames(Y_subset)[i], "~", colnames(Y_subset)[j]))
      model <- pgls(formula, data = compdata)
      summary_model <- summary(model)
      pgls_p_values[i, j] <- summary_model$coefficients[2, 4]  # p-value for the predictor
      pgls_adj_R2_values[i, j] <- summary_model$adj.r.squared  # Adjusted R^2 value
    } else {
      pgls_p_values[i, j] <- NA
      pgls_adj_R2_values[i, j] <- NA
    }
  }
}

# Combine p-values and adjusted R^2 into a single matrix for heatmap
combined_matrix <- matrix(paste0(
  "p: ", format(pgls_p_values, digits = 2), "\n",
  "adj. R^2: ", format(pgls_adj_R2_values, digits = 2)
), nrow=nrow(pgls_p_values), ncol=ncol(pgls_p_values))
rownames(combined_matrix) <- colnames(Y_subset)
colnames(combined_matrix) <- colnames(Y_subset)

# Define color palette for heatmap
col_palette <- colorRampPalette(c("white", "red"))(200)

# Create a heatmap for p-values and adjusted R^2 values
corrplot(
  pgls_adj_R2_values, method="color", type="upper",
  addCoef.col = "black", # Add adjusted R^2 values
  tl.col="black", tl.srt=45, # Text label color and rotation
  title = "Heatmap of Adjusted R^2 Values from PGLS Models",
  col = col_palette, is.corr = FALSE, na.label = "NA"
)

# Create a heatmap for p-values
corrplot(
  pgls_p_values, method="color", type="upper",
  addCoef.col = "black", # Add p-values
  tl.col="black", tl.srt=45, # Text label color and rotation
  title = "Heatmap of p-values from PGLS Models",
  col = col_palette, is.corr = FALSE, na.label = "NA"
)

# Combine p-values and adjusted R^2 into a single matrix for visualization
combined_matrix <- matrix(paste0(
  "p: ", format(pgls_p_values, digits = 2), "\n",
  "adj. R^2: ", format(pgls_adj_R2_values, digits = 2)
), nrow=nrow(pgls_p_values), ncol=ncol(pgls_p_values))
rownames(combined_matrix) <- colnames(Y_subset)
colnames(combined_matrix) <- colnames(Y_subset)

# Display the combined matrix
print(combined_matrix)

# Create a heatmap for the combined matrix
heatmap(
  as.matrix(combined_matrix), Rowv = NA, Colv = NA, 
  col = col_palette, scale = "none", margins = c(5, 10)
)

# Load necessary libraries
library(caper)
library(ape)
library(ggplot2)
library(gridExtra)

# Define PGLS models
pgls_bamm <- pgls(DR ~ BAMM, data = compdata, lambda = "ML")
pgls_clads <- pgls(DR ~ ClaDs, data = compdata, lambda = "ML")
pgls_misse <- pgls(DR ~ MISSE, data = compdata)
pgls_bamm_clads <- pgls(BAMM ~ ClaDs, data = compdata, lambda = "ML")
pgls_bamm_misse <- pgls(BAMM ~ MISSE, data = compdata, lambda = "ML")
pgls_clads_misse <- pgls(ClaDs ~ MISSE, data = compdata, lambda ="ML")

# Extract residuals and fitted values
extract_residuals_fitted <- function(model) {
  residuals <- residuals(model)
  fitted <- fitted(model)
  list(residuals = residuals, fitted = fitted)
}

res_fit_bamm <- extract_residuals_fitted(pgls_bamm)
res_fit_clads <- extract_residuals_fitted(pgls_clads)
res_fit_misse <- extract_residuals_fitted(pgls_misse)
res_fit_bamm_clads <- extract_residuals_fitted(pgls_bamm_clads)
res_fit_bamm_misse <- extract_residuals_fitted(pgls_bamm_misse)
res_fit_clads_misse <- extract_residuals_fitted(pgls_clads_misse)

# Create residuals vs fitted values plots and QQ plots using ggplot2
create_plots <- function(fitted, residuals, model_name) {
  p1 <- ggplot(data.frame(Fitted = fitted, Residuals = residuals), aes(Fitted, Residuals)) +
    geom_point() +
    geom_smooth(method = "loess", se = FALSE, color = "red") +
    ggtitle(paste(model_name, ": Residuals vs Fitted"))
  
  p2 <- ggplot(data.frame(Sample = residuals), aes(sample = Sample)) +
    stat_qq() +
    stat_qq_line() +
    ggtitle(paste(model_name, ": Normal Q-Q Plot"))
  
  list(p1 = p1, p2 = p2)
}

plots_bamm <- create_plots(res_fit_bamm$fitted, res_fit_bamm$residuals, "DR ~ BAMM")
plots_clads <- create_plots(res_fit_clads$fitted, res_fit_clads$residuals, "DR ~ ClaDS")
plots_misse <- create_plots(res_fit_misse$fitted, res_fit_misse$residuals, "DR ~ MISSE")
plots_bamm_clads <- create_plots(res_fit_bamm_clads$fitted, res_fit_bamm_clads$residuals, "BAMM ~ ClaDS")
plots_bamm_misse <- create_plots(res_fit_bamm_misse$fitted, res_fit_bamm_misse$residuals, "BAMM ~ MISSE")
plots_clads_misse <- create_plots(res_fit_clads_misse$fitted, res_fit_clads_misse$residuals, "ClaDS ~ MISSE")

# Arrange the plots in a grid
quartz()
grid.arrange(plots_bamm$p1, plots_bamm$p2,
             plots_clads$p1, plots_clads$p2,
             plots_misse$p1, plots_misse$p2,
             plots_bamm_clads$p1, plots_bamm_clads$p2,
             plots_bamm_misse$p1, plots_bamm_misse$p2,
             plots_clads_misse$p1, plots_clads_misse$p2,
             ncol = 4)

install.packages("mice")
Yinstall.packages("naniar")
library(mice)
library(naniar)

data <- syg[, c("Length","Clutch_size","Egg_size","Brood_size")]

# Summarize the pattern of missing data
md_pattern <- md.pattern(data)
quartz()
print(md_pattern)

# Visualize the missing data pattern
library(ggplot2)
gg_miss_upset(data)

# Perform logistic regression test for MCAR
mcar_test <- mice::mcar(data)
print(mcar_test)
