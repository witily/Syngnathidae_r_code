library(scales)
library(ggplot2)
library(tidyr)
library(ggtree)
library(ggpubr)
library(reshape2)
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
load("/Users/jeremyhappe/Desktop/Condorcet/TFE/syg.RData")
tree <- read.tree(file = "/Users/jeremyhappe/Desktop/Condorcet/TFE/BAMM_new_tree/full_tree.tre")
#############
Y <- syg
Y <- na.omit(Y[, c("ID", "Species","Strips")])
common_taxa <- intersect(tree$tip.label, Y$Species)
if (length(common_taxa) == 0) {
  stop("No common taxa found between the tree and Excel file.")}
tree_pruned <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% common_taxa])

Y_filtered <- Y %>% filter(Species %in% common_taxa)

# Ensure the order of species in Y_filtered matches the order of tips in the pruned tree
Y_filtered <- Y_filtered %>% arrange(match(Species, tree_pruned$tip.label))

# Prepare data for plotting
trait_data <- data.frame(Species = Y_filtered$Species, Strips = as.factor(Y_filtered$Strips))

# Convert tree to ggtree object
ggtree_obj <- ggtree(tree_pruned) %<+% trait_data + 
  geom_tiplab(aes(label = Species, color = Strips), show.legend = FALSE) +
  scale_color_manual(values = c("0" = "red", "1" = "green")) +
  theme_tree() +
  labs(color = "Strips Presence (1) / Absence (0)") +
  theme(legend.position = "right")

# Print the plot
print(ggtree_obj)




#rename 
Y$BAMM <- Y$speciation_rate

#colour designation
families <- c("Aulostomidae", "Callionymidae", "Centriscidae", "Dactylopteridae", "Fistulariidae", "Mullidae", "Pegasidae", "Solenostomidae", "Syngnathidae")
colors <- c("#F8766D", "#D39200", "#93AA00", "#00BA38" ,"#00C19F", "#00B9E3", "#619CFF" ,"#DB72FB" ,"#FF61C3")
family_colors <- setNames(colors, families)

## visualize the tree 

p <- ggtree(tree,layout = 'rect', right =TRUE) 
p <- p %<+% Y + geom_tippoint(aes(color=Family), size = 1, show.legend = FALSE) + geom_tiplab(size = 2.6,as_ylab = TRUE,hjust = -0.8 )
p

#Plot with  categories
library(tidyverse)
str(Y)
spe <- pivot_longer(Y, c(DR, speciation_rate),values_to = "bob", names_to = "Speciation_rate")

p <- facet_plot(p, data = spe, geom=geom_bar,
                mapping=aes(x=bob,fill= Speciation_rate),
                width=0.7,
                stat="identity",
                orientation="y",
                panel = "Speciation rate") + scale_fill_manual(values=c("DR" = "deepskyblue3", "speciation_rate"="darkolivegreen2")) + xlim_tree(2.1)+
  theme_tree2() + 
  theme(
    legend.position.inside = c(.85, .65),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12)  
  )

p +geom_facet(panel = "speciation rate", data = Y, geom = geom_col, 
             aes(x = dummy_bar_value, color = Family, 
                 fill = Family), orientation = 'y', width = .6) +
  theme_tree2(legend.position=c(.05, .85))

# Reshape from wide to long format
Y$Sexual_ornements <-as.factor(Y$Sexual_ornements)
Y$Courtship <-as.factor(Y$Courtship)
Y$Color_dimorphism<-as.factor(Y$Color_dimorphism)
Y$Colour_dimorphism <-Y$Color_dimorphism

#PLOT Cumulative bars

#sexual traits
long_data_color <- Y %>%
  pivot_longer(cols = c("Courtship", "Monogamous","Sexual_dimorphism","Sexual_ornements"
  ),
  names_to = "Trait", values_to = "Value")

#colouration traits
long_data_color <- Y %>%
  pivot_longer(cols = c("Colour_dimorphism", "Red",
                        "Conspicuous", "Strips", "Dots",
                        "Mimicry", "Camouflage", "Polymorphic"
                         ),
               names_to = "Trait", values_to = "Value") 

percent_data <- long_data_color %>%
  group_by(Family, Trait, Value) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Family, Trait) %>%
  mutate(Percentage = Count / sum(Count) * 100)

ggplot(percent_data, aes(x = Family, y = Percentage, fill = as.factor(Value))) +
  geom_bar(stat = "identity", position = "stack",size = 0.4, color = "black") +
  geom_text(aes(label = ifelse(is.na(Value), "", sprintf("%.1f%%", Percentage))), 
            position = position_stack(vjust = 0.5), 
            color = "white",
            size = 1.9) +  
  facet_wrap(~ Trait, scales = "free_x", nrow = 2, ) +
  labs(x = "Family", y = "Percentage (%)", fill = "Binary Value", title = "Cumulative Distribution of Traits by Family") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.margin = unit(c(1,1,1,1), "cm"),
        text = element_text(size = 8)) +  # General text size for other plot elements can also be adjusted here
  scale_fill_manual(values = c("0" = "#619CFF", "1" = "orange"))



# Assuming your data frame's species column is named 'Species'
species_count <-Y %>%
  count(Family) %>%
  arrange(desc(n))  # This arranges the species by the number of individuals, descending
library(forcats)  # Pour fct_rev

species_count$Family <- fct_rev(species_count$Family)

#NUMBER OF SPECIES

ggplot(species_count, aes(x = Family, y = n, fill = Family)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Rendre les barres horizontales
  geom_text(aes(label = n), hjust = -0.5, color = "black", size = 4.5) +  # Ajuster les paramètres de texte
  labs(x = "Species", y = "Family", title = "Family") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12, angle = 0),  # Agrandir et horizontaliser les noms des familles
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.margin = unit(c(1, 1, 1, 1), "cm")
       ) +  # Agrandir la légende
  scale_fill_manual(values = family_colors)  # Utiliser l'échelle de couleur personnalisée


#barplot for continous trait

Y$Brood_size <- exp(Y$Brood_size)
Y$Egg_size <- exp(Y$Egg_size)
Y$Clutch_size <- exp(Y$Clutch_size)

# Named vector for custom facet labels
facet_labels <- c(Brood_size = "Brood Size", Clutch_size = "Clutch Size", Egg_size = "Egg Size (mm)")

# Plot
p <- ggplot(Y_long, aes(x = Family, y = Size, fill = Trait)) +
  geom_boxplot(alpha = 0.8) +  # Adjust the transparency of the boxplot
  facet_wrap(~Trait, scales = "free", labeller = labeller(Trait = facet_labels)) +  
  labs(title = "Box Plots of Egg Size, Brood Size, and Clutch Size by Family", x = "Family", y = NULL) +
  scale_y_log10() +  # Apply log scale to the y-axis
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis text for better readability
        legend.title = element_blank(), # Remove the legend title if it's redundant
        strip.background = element_rect(fill = "white", color = "black", size = 1), # Add frames around facets
        strip.text = element_text(face = "bold")) # Make the facet labels bold


# NOT USED

Y$Length <- exp(Y$Length)

p <- ggtree(tree,right = TRUE) 
# Add the bar plot next to the tree
p <- facet_plot(p, panel = "Strips", data = Y, geom = geom_bar, 
                mapping = aes(x = Length, fill = as.factor(ID)), 
                width = 0.7, stat = "identity", orientation = "y") +
  scale_fill_manual(values = rep("deepskyblue3", length(unique(Y$Species)))) +  # Use one color
  theme_tree2() +
  theme(
    legend.position = "none"  # Remove legend since it's not needed
  )

# Print the plot
print(p)

ggsave(filename = "tree_plot_l.png", plot = p, width = 10, height = 14, dpi = 300)  # Adjust width, height, and dpi as needed


library(DescTools)
library(Hmisc)

p <- ggtree(tree, layout = 'rect', right = TRUE) 
# Add the bar plot next to the tree with different colors for each family
p <- facet_plot(p, panel = "Length", data = Y, geom = geom_bar, 
                mapping = aes(x = Length, fill = Family),  # Map fill to Family
                width = 0.7, stat = "identity", orientation = "y") +
  scale_x_log10()+
  scale_fill_manual(values = family_colors) +  # Use custom color scale
  theme_tree2() +
  theme(
    legend.position = "none"  # Remove legend if not needed
  )
p
p <- ggtree(tree, right = TRUE)
library(skimr)
skim(Y)

Y$Clutch_size <- log(Y$Clutch_size)
Y$Brood_size <- log(Y$Brood_size)
Y$Egg_size <- log(Y$Egg_size)
# Tracer la relation entre Length et Clutch_size
p1 <- ggplot(Y, aes(x = Length, y = Clutch_size)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed") +
  labs(title = "Relation between Length and Clutch_size",
       x = "Length",
       y = "Clutch_size")

# Tracer la relation entre Length et Egg_size
p2 <- ggplot(Y, aes(x = Length, y = Egg_size)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed") +
  labs(title = "Relation between Length and Egg_size",
       x = "Length",
       y = "Egg_size")


p1 <- ggplot(Y, aes(x = Clutch_size, y = Brood_size)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), method = "pearson", label.x = 3) + # Pearson correlation and R2
  stat_cor(aes(label = paste("Spearman: ", ..r.label.., sep = "")), method = "spearman", label.x = 3, label.y = max(Y$Brood_size) * 0.9) + # Spearman correlation
  labs(title = "Relation between Length and Brood_size",
       x = "Length",
       y = "Brood_size")

# Tracer la relation entre Length et Egg_size
p2 <- ggplot(Y, aes(x = Clutch_size, y = Egg_size)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), method = "pearson", label.x = 3) + # Pearson correlation and R2
  stat_cor(aes(label = paste("Spearman: ", ..r.label.., sep = "")), method = "spearman", label.x = 3, label.y = max(Y$Egg_size) * 0.9) + # Spearman correlation
  labs(title = "Relation between Clutch_size and Egg_size",
       x = "Clutch_size",
       y = "Egg_size")

# Afficher les graphiques
print(p1)
print(p2)

# Calculer les corrélations de Spearman et Pearson
cor_length_clutch_spearman <- cor(Y$Length, Y$Clutch_size, method = "spearman")
cor_length_clutch_pearson <- cor(Y$Length, Y$Clutch_size, method = "pearson")
cor_length_egg_spearman <- cor(Y$Length, Y$Egg_size, method = "spearman")
cor_length_egg_pearson <- cor(Y$Length, Y$Egg_size, method = "pearson")

# Afficher les résultats de corrélation
cat("Spearman correlation between Length and Clutch_size:", cor_length_clutch_spearman, "\n")
cat("Pearson correlation between Length and Clutch_size:", cor_length_clutch_pearson, "\n")
cat("Spearman correlation between Length and Egg_size:", cor_length_egg_spearman, "\n")
cat("Pearson correlation between Length and Egg_size:", cor_length_egg_pearson, "\n")



#Body size TUSKEY test
Y$Length<-exp(Y$Length)
Y$Length <- log(Y$Length)
ggplot(Y, aes(x = Family, y = Length)) +
  geom_boxplot(aes(fill=Family)) +
  labs(title = "Boxplot of Body size by Family",
       x = "Family",
       y = "Body size (cm)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

anova_result <- aov(Length ~ Family, data = Y)
summary(anova_result)


tukey_result <- TukeyHSD(anova_result)
tukey_result_df <- as.data.frame(tukey_result$Family)
tukey_result_df$comparison <- rownames(tukey_result_df)

heatmap_data <- tukey_result_df %>%
  mutate(significant = ifelse(`p adj` < 0.05, "Significant", "Not Significant")) %>%
  separate(comparison, into = c("group1", "group2"), sep = "-")

# Reshape data for heatmap
heatmap_matrix <- acast(heatmap_data, group1 ~ group2, value.var = "significant", fill = "Not Compared")

# Create the heatmap
ggplot(melt(heatmap_matrix), aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("Not Significant" = "blue", "Significant" = "red", "Not Compared" = "grey")) +
  labs(title = "Significance Table of Family Size Differences based on a TuskeyHSD test",
       x = "Family 1",
       y = "Family 2",
       fill = "Significance") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )+
  coord_flip()

