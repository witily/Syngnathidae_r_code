library(hisse)

load("/Users/jeremyhappe/Desktop/Condorcet/TFE/syg.RData")
tree <- read.tree(file = "/Users/jeremyhappe/Desktop/Condorcet/TFE/BAMM_new_tree/full_tree.tre")

f <- ((3+33+7+6+4+51+2+4+159)/(3+201+13+7+4+104+8+6+325)) #Conservative
f <- ((3+33+7+6+4+51+2+4+159)/(3+201+13+7+4+104+8+6+338)) #"new species" inculded in fractions

#can't run it 
choice <-MiSSEGreedy(tree, f=f,n.cores= 4, save.file = "/Users/jeremyhappe/Desktop/Condorcet/TFE/")

####Single rate model
turnover <- c(1)
eps <- c(1)
one.rate_f <- MiSSE(tree, f=f, turnover=turnover, eps=eps)
#### Double rate model
turnover <- c(1,2)
eps <- c(1,1)
two.rate_f <- MiSSE(tree, f=f, turnover=turnover, eps=eps)
#rate classes A:C
turnover <- c(1,2,3)
eps <- c(1,1,1)
three.rate_f <- MiSSE(tree, f=f, turnover=turnover, eps=eps)
#rate classes A:D
turnover <- c(1,2,3,4)
eps <- c(1,1,1,1)
four.rate_f <- MiSSE(tree, f=f, turnover=turnover, eps=eps)
#rate classes A:E
turnover <- c(1,2,3,4,5)
eps <- c(1,1,1,1,1)
five.rate_f <- MiSSE(tree, f=f, turnover=turnover, eps=eps)
#rate classes A:F
turnover <- c(1,2,3,4,5,6)
eps <- c(1,1,1,1,1,1)
six.rate <- MiSSE(tree, f=1, turnover=turnover, eps=eps)

one.rate_f$AIC
two.rate_f$AIC
three.rate_f$AIC
four.rate_f$AIC
five.rate_f$AIC


save(two.rate_f, one.rate_f, three.rate_f, four.rate_f,five.rate_f, file = "/Users/jeremyhappe/Desktop/Condorcet/TFE/saved_rates_f.RData")
# load("saved_rates.RData")

load("/Users/jeremyhappe/Desktop/Condorcet/TFE/saved_rates_f.RData")

three.rate.recon <- MarginReconMiSSE(phy=tree, f=1, hidden.states=2,
                                     pars=three.rate_f$solution, n.cores=3, AIC=two.rate_f$AIC)

three.rate.recon <- MarginReconMiSSE(phy=tree, f=1, hidden.states=3,
                                   pars=three.rate_f$solution, n.cores=3, AIC=three.rate_f$AIC)

save(two.rate.recon, three.rate.recon, file = "/Users/jeremyhappe/Desktop/Condorcet/TFE/models_MISSE.RData")
load("/Users/jeremyhappe/Desktop/Condorcet/TFE/models_MISSE.RData")

plot.misse.states(two.rate.recon, rate.param="speciation", show.tip.label=TRUE, type="phylogram",
                  fsize=.25, legend= "tips")
plot.misse.states(three.rate.recon, rate.param="speciation", show.tip.label=TRUE, type="phylogram",
                  fsize=.25, legend= "tips")

speciation_rates <- GetModelAveRates(misse.results.list, type = c("tips"))

save(speciation_rates, file = "/Users/jeremyhappe/Desktop/Condorcet/TFE/MISSE_speciation.RData")

misse.results.list = list()
misse.results.list[[1]] = two.rate.recon
misse.results.list[[2]] = three.rate.recon

plot.misse.states(misse.results.list, rate.param="speciation", show.tip.label=TRUE, type="phylogram",
                  fsize=.25, legend="tips")
?plot.misse.states

print(three.rate.recon$rates.mat)

## Vizualise 

aic_values <- c(
  "1.rate" = 2036.901,
  "2.rate" = 2010.608,
  "3.rate" = 2002.738,
  "4.rate" = 2004.197,
  "5.rate" = 2010.509
)


# Create the line plot
aic_df <- data.frame(
  rate_f = factor(names(aic_values)),
  AIC = as.numeric(aic_values)
)

# Create the line plot
library(ggplot2)

ggplot(aic_df, aes(x = rate_f, y = AIC, group = 1)) +
  geom_line() +
  geom_point() +
  labs(title = "AIC Values for Different Rate Models",
       x = "Models",
       y = "AIC") +
  theme_minimal(base_size = 15) +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  coord_cartesian(ylim = c(2000, 2040))




