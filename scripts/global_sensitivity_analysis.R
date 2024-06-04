# global sensitivity analysis using random forest

# read in packages
library(ggplot2)
library(deSolve)
library(rootSolve)
library(reshape2)
library(dplyr)
library(ggpol) 
library(MetBrewer) 
library(ggpubr)
library(randomForest)

# read in parameter simulations
param_viable <- read.csv("output/param_simulations_stable_10000.csv") %>% head(1000)

# create a new column of the mean of all compartment values
param_viable <- param_viable %>%
mutate(model = (Pt+Nt+Lt+Ha+Pa+Na)/6)

# define a vector of state variables to loop through
state_variables <- c("Pt", "Nt", "Lt", "Ha", "Pa", "Na", "pt_prod", "ha_prod", "pa_prod", "model")

# allocate an empty dataframe to contain the results of each sensitivity analysis 
sensitivity_results <- NULL

# run sensitivity analysis on all state variables
for (i in 1:length(state_variables)) {

# isolate columns to work with
sensitivity_columns <- param_viable %>% select(state_variables[i], αt:la)

# rename the state variable column to a neutral name
colnames(sensitivity_columns)[1] <- "compartment"

# run the randomForest model to estimate importance of each parameter
sa_temp <- randomForest(x=(sensitivity_columns %>% select(-compartment)),
 y=sensitivity_columns$compartment, 
 proximity = TRUE, 
 importance = TRUE)

# extract the importance values of each parameter
importance_temp <- data.frame(sa_temp[["importance"]]) %>% 
 rownames_to_column(., var = "parameter") %>% 
 mutate(compartment = state_variables[i])

# now calculate the ranked importance value as the ratio of each parameter's 
# importance values to the sum of all parameters' importance values
ranked_importance_temp <- importance_temp %>% 
 mutate(ranked_importance = IncNodePurity/sum(importance_temp$IncNodePurity))

# store the results from each loop run in the sensitivity_results dataframe
sensitivity_results <- bind_rows(sensitivity_results, ranked_importance_temp)

}

# plot the results
gsa_plot <- ggplot(sensitivity_results, aes(x=parameter, y=ranked_importance, fill=parameter)) +
geom_col(alpha=0.9) + 
theme_classic() + 
xlab("parameter") + 
ylab("relative importance") + 
facet_wrap(vars(compartment), ncol = 4) + 
scale_fill_viridis_d() +
  guides(fill=FALSE)

ggsave(file = "output/gsa_plot.svg", plot=gsa_plot, width=14, height=6)
write.csv(sensitivity_results, "output/gsa_results.csv", row.names=FALSE)

 
