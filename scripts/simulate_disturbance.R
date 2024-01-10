# read in packages
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(tictoc)
library(lhs)
library(data.table)

# generate a dataset with numbers ranging from 0-10
param_estimates1 <- NULL
for(i in seq(1,4000,1)){
  
  # sample 100 times for 15 parameters
  loop_1 <- randomLHS(100, 17)
  
  # bind each iteration on top of the other 
  # each column is a list of 10000 samples with equal distribution in 100 bins
  # multiplied by 10 to get range of data from 0 to 10
  param_estimates1 <- rbind(param_estimates1, data.frame(Result=loop_1*10))
}

# generate a dataset with numbers ranging from 0-1
param_estimates2 <- NULL
for(i in seq(1,4000,1)){
  ## sample 100 times for 6 parameters
  loop_2 <- randomLHS(100, 7)
  
  # bind each iteration on top of the other 
  # each column is a list of 10000 samples with equal distribution in 100 bins
  param_estimates2 <- rbind(param_estimates2, data.frame(Result=loop_2))
}

# export files (for records)
write.csv(param_estimates1, "output/param_estimates1.csv", row.names=FALSE)
write.csv(param_estimates2, "output/param_estimates2.csv", row.names=FALSE)

tic("generate parameter combinations")
# Next, make a dataframe for 1000 parameter combinations
columns <- c("αt", 
             "αa",
             "βa",
             "e",
             "ρ",
             "δ",
             "γ",
             "ϵ",
             "ηa",
             "θt",
             "θa",
             "λt",
             "λa",
             "μt",
             "μa",
             "τa",
             "ψt",
             "ψa",
             "lt",
             "la",
             "Pt",
             "Nt",
             "Lt",
             "Ha",
             "Pa",
             "Na",
             "ha_prod",
             "pt_prod",
             "pa_prod",
             "pt_flux",
             "nt_flux",
             "na_flux")

param_simulations <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(param_simulations) = columns

# generate a dataframe of 10000 parameter simulations
for (i in seq(400000)){ 
  
  # select rows to draw from
  params_1 <- param_estimates1[i,]
  params_2 <- param_estimates2[i,]
  
  # producer uptake rates
  αt <- params_1[[1]] 
  αa <- params_1[[2]] 
  
  # herbivore uptake rates
  βa <- params_1[[3]] 
  e <- params_2[[5]] 
  ρ <- params_2[[6]] 
  
  # producer recycling rates 
  μt <- params_2[[1]] 
  μa <- params_2[[2]] 
  
  # herbivore recycling rates 
  ηa <- params_2[[4]] 
  
  # producer emigration rates (leaving the system)
  θt <- params_1[[6]]
  θa <- params_1[[7]]
  
  # herbivore emigration rates (leaving the system)
  τa <- params_1[[8]]
  
  # rate of movement between Pt and aquatic ecosystem
  ϵ <- params_1[[9]]
  γ <- params_1[[10]]
  
  # proportion of Pt partitioned to Ha rather than Na
  δ <- params_1[[17]]
  
  # flow rate between inorganic nitrogen pools
  ψt <- params_1[[11]] 
  ψa <- params_1[[12]] 
  
  # input to inorganic nitrogen pools
  λt <- params_1[[13]]
  λa <- params_1[[14]]
  
  # leaching rate from nitrogen pools
  lt <- params_1[[15]]
  la <- params_1[[16]]
  
  # equilibrium values and ecosystem properties
  # equation 10
  Pt <- -0.5*(-2*e*la*lt*αa*βa*δ*(ϵ**2) - 4*e*la*lt*αa*βa*δ*ϵ*θt -  2*e*la*lt*αa*βa*δ*(θt**2) + 2*e*la*αa*αt*βa*δ*ϵ*λt +  
                2*e*la*αa*αt*βa*δ*θt*λt + 2*e*la*lt*αa*βa*δ*ϵ*θt*μt +  2*e*la*lt*αa*βa*δ*(θt**2)*μt - 2*e*la*αa*αt*βa*δ*θt*λt*μt -  
                2*lt*(αa**2)*δ*(ϵ**2)*τa + 2*e*lt*(αa**2)*δ*(ϵ**2)*ηa*τa -  4*lt*(αa**2)*δ*ϵ*θt*τa + 4*e*lt*(αa**2)*δ*ϵ*ηa*θt*τa -  
                2*lt*(αa**2)*δ*(θt**2)*τa + 2*e*lt*(αa**2)*δ*ηa*(θt**2)*τa +  2*(αa**2)*αt*δ*ϵ*λt*τa - 2*e*(αa**2)*αt*δ*ϵ*ηa*λt*τa +  
                2*(αa**2)*αt*δ*θt*λt*τa - 2*e*(αa**2)*αt*δ*ηa*θt*λt*τa +  2*lt*(αa**2)*δ*ϵ*θt*μt*τa - 2*e*lt*(αa**2)*δ*ϵ*ηa*θt*μt*τa +  
                2*lt*(αa**2)*δ*(θt**2)*μt*τa -  2*e*lt*(αa**2)*δ*ηa*(θt**2)*μt*τa - 2*(αa**2)*αt*δ*θt*λt*μt*τa +  2*e*(αa**2)*αt*δ*ηa*θt*λt*μt*τa - 
                e*la*αt*(βa**2)*γ*ϵ*ψa -  2*e*lt*αa*βa*δ*(ϵ**2)*ψa + e*la*αt*βa*δ*ϵ*θa*ψa -  e*la*αt*(βa**2)*γ*θt*ψa - 4*e*lt*αa*βa*δ*ϵ*θt*ψa +  
                e*la*αt*βa*δ*θa*θt*ψa - 2*e*lt*αa*βa*δ*(θt**2)*ψa +  e*αa*αt*βa*δ*ϵ*λa*ψa + e*αa*αt*βa*δ*θt*λa*ψa + 2*e*αa*αt*βa*δ*ϵ*λt*ψa +  
                2*e*αa*αt*βa*δ*θt*λt*ψa + e*la*αt*(βa**2)*γ*θt*μt*ψa +  2*e*lt*αa*βa*δ*ϵ*θt*μt*ψa - e*la*αt*βa*δ*θa*θt*μt*ψa +  
                2*e*lt*αa*βa*δ*(θt**2)*μt*ψa - e*αa*αt*βa*δ*θt*λa*μt*ψa -  2*e*αa*αt*βa*δ*θt*λt*μt*ψa + lt*αa*βa*δ*(ϵ**2)*ρ*ψa +  
                lt*αa*βa*δ*ϵ*θt*ρ*ψa - αa*αt*βa*δ*ϵ*λt*ρ*ψa - αa*αt*βa*γ*ϵ*τa*ψa +  e*αa*αt*βa*γ*ϵ*ηa*τa*ψa + αa*αt*δ*ϵ*θa*τa*ψa - 
                2*e*αa*αt*δ*ϵ*ηa*θa*τa*ψa -  αa*αt*βa*γ*θt*τa*ψa + e*αa*αt*βa*γ*ηa*θt*τa*ψa + αa*αt*δ*θa*θt*τa*ψa -  
                2*e*αa*αt*δ*ηa*θa*θt*τa*ψa + αa*αt*δ*ϵ*θa*μa*τa*ψa + αa*αt*δ*θa*θt*μa*τa*ψa +  αa*αt*βa*γ*θt*μt*τa*ψa - 
                e*αa*αt*βa*γ*ηa*θt*μt*τa*ψa -  αa*αt*δ*θa*θt*μt*τa*ψa + 2*e*αa*αt*δ*ηa*θa*θt*μt*τa*ψa -  αa*αt*δ*θa*θt*μa*μt*τa*ψa + 
                e*αt*βa*δ*ϵ*θa*(ψa**2) -  e*αt*(βa**2)*γ*θt*(ψa**2) + e*αt*βa*δ*θa*θt*(ψa**2) +  e*αt*(βa**2)*γ*θt*μt*(ψa**2) - 
                e*αt*βa*δ*θa*θt*μt*(ψa**2) -  αt*βa*δ*ϵ*θa*μa*ρ*(ψa**2) - 2*e*la*αa*βa*δ*(ϵ**2)*ψt -  4*e*la*αa*βa*δ*ϵ*θt*ψt - 
                2*e*la*αa*βa*δ*(θt**2)*ψt +  2*e*la*αa*βa*δ*ϵ*θt*μt*ψt + 2*e*la*αa*βa*δ*(θt**2)*μt*ψt -  2*(αa**2)*δ*(ϵ**2)*τa*ψt + 
                2*e*(αa**2)*δ*(ϵ**2)*ηa*τa*ψt -  4*(αa**2)*δ*ϵ*θt*τa*ψt + 4*e*(αa**2)*δ*ϵ*ηa*θt*τa*ψt -  2*(αa**2)*δ*(θt**2)*τa*ψt + 
                2*e*(αa**2)*δ*ηa*(θt**2)*τa*ψt +  2*(αa**2)*δ*ϵ*θt*μt*τa*ψt - 2*e*(αa**2)*δ*ϵ*ηa*θt*μt*τa*ψt +  2*(αa**2)*δ*(θt**2)*μt*τa*ψt -  
                2*e*(αa**2)*δ*ηa*(θt**2)*μt*τa*ψt - e*αa*βa*δ*(ϵ**2)*ψa*ψt -  2*e*αa*βa*δ*ϵ*θt*ψa*ψt - e*αa*βa*δ*(θt**2)*ψa*ψt +  
                e*αa*βa*δ*ϵ*θt*μt*ψa*ψt + e*αa*βa*δ*(θt**2)*μt*ψa*ψt +  αa*βa*δ*(ϵ**2)*ρ*ψa*ψt + αa*βa*δ*ϵ*θt*ρ*ψa*ψt +  
                ψa*sqrt(4*αa*αt*δ*(αa*(ϵ + θt - θt*μt)*τa - βa*ϵ*ρ*ψa + e*(ϵ + θt - θt*μt)*(la*βa - αa*ηa*τa + βa*ψa))* 
                          (-(e*(lt*(βa**2)*γ*ϵ*(ϵ + θt) -  αt*((βa**2)*γ*(ϵ*(λa + λt) - θt*λa*(-1 + μt)) +  δ*ηa*(θa**2)*
                                                                 (ϵ + θt - θt*μt)*τa -  βa*θa*(ϵ + θt - θt*μt)*(δ*λa + γ*ηa*τa)) + 
                                  βa*(ϵ + θt)*(βa*γ*θt*(-1 + μt) + δ*θa*(ϵ + θt - θt*μt))*ψt)) + 
                             θa*μa*(lt*βa*δ*ϵ*(ϵ + θt)*ρ -αt*(βa*δ*ϵ*λt*ρ - βa*γ*(ϵ + θt - θt*μt)*τa +δ*θa*(ϵ + θt - θt*μt)*τa) + 
                                      βa*δ*ϵ*(ϵ + θt)*ρ*ψt)) + (lt*αa*βa*δ*ϵ*(ϵ + θt)*ρ - αa*αt*βa*δ*ϵ*λt*ρ + αa*αt*βa*γ*ϵ*τa -αa*αt*δ*ϵ*θa*τa + 
                                                                  αa*αt*βa*γ*θt*τa - αa*αt*δ*θa*θt*τa -αa*αt*δ*ϵ*θa*μa*τa - αa*αt*δ*θa*θt*μa*τa - 
                                                                  αa*αt*βa*γ*θt*μt*τa +αa*αt*δ*θa*θt*μt*τa + αa*αt*δ*θa*θt*μa*μt*τa + αt*βa*δ*ϵ*θa*μa*ρ*ψa +
                                                                  αa*βa*δ*(ϵ**2)*ρ*ψt + αa*βa*δ*ϵ*θt*ρ*ψt +
                                                                  e*(la*αt*βa*(βa*γ - δ*θa)*(ϵ + θt - θt*μt) -   αt*βa*(βa*γ*θt*(-1 + μt) + δ*θa*(ϵ + θt - θt*μt))*ψa -   
                                                                       αa*(ϵ + θt - θt*μt)*   (αt*(βa*δ*λa + βa*γ*ηa*τa - 2*δ*ηa*θa*τa) + βa*δ*(ϵ + θt)*ψt))**2)))/ 
    (αa*αt*δ*(ϵ + θt - θt*μt)*(-(αa*(ϵ + θt - θt*μt)*τa) + βa*ϵ*ρ*ψa -e*(ϵ + θt - θt*μt)*(la*βa - αa*ηa*τa + βa*ψa)))
  Nt <- (ϵ + θt)/αt
  Lt <- -0.5*(-(e*la*αt*(βa**2)*γ*ϵ) + e*la*αt*βa*δ*ϵ*θa - e*la*αt*(βa**2)*γ*θt + e*la*αt*βa*δ*θa*θt - e*αa*αt*βa*δ*ϵ*λa - e*αa*αt*βa*δ*θt*λa + e*la*αt*(βa**2)*γ*θt*μt - e*la*αt*βa*δ*θa*θt*μt + e*αa*αt*βa*δ*θt*λa*μt + lt*αa*βa*δ*(ϵ**2)*ρ + lt*αa*βa*δ*ϵ*θt*ρ - αa*αt*βa*δ*ϵ*λt*ρ - αa*αt*βa*γ*ϵ*τa + e*αa*αt*βa*γ*ϵ*ηa*τa + αa*αt*δ*ϵ*θa*τa - αa*αt*βa*γ*θt*τa + e*αa*αt*βa*γ*ηa*θt*τa + αa*αt*δ*θa*θt*τa - αa*αt*δ*ϵ*θa*μa*τa - αa*αt*δ*θa*θt*μa*τa + αa*αt*βa*γ*θt*μt*τa - e*αa*αt*βa*γ*ηa*θt*μt*τa - αa*αt*δ*θa*θt*μt*τa + αa*αt*δ*θa*θt*μa*μt*τa + e*αt*βa*δ*ϵ*θa*ψa - e*αt*(βa**2)*γ*θt*ψa + e*αt*βa*δ*θa*θt*ψa + e*αt*(βa**2)*γ*θt*μt*ψa - e*αt*βa*δ*θa*θt*μt*ψa - αt*βa*δ*ϵ*θa*μa*ρ*ψa - e*αa*βa*δ*(ϵ**2)*ψt - 2*e*αa*βa*δ*ϵ*θt*ψt - e*αa*βa*δ*(θt**2)*ψt + e*αa*βa*δ*ϵ*θt*μt*ψt + e*αa*βa*δ*(θt**2)*μt*ψt + αa*βa*δ*(ϵ**2)*ρ*ψt + αa*βa*δ*ϵ*θt*ρ*ψt + sqrt(4*αa*αt*δ*(αa*(ϵ + θt - θt*μt)*τa - βa*ϵ*ρ*ψa + e*(ϵ + θt - θt*μt)*(la*βa - αa*ηa*τa + βa*ψa))*(-(e*(lt*(βa**2)*γ*ϵ*(ϵ + θt) - αt*((βa**2)*γ*(ϵ*(λa + λt) - θt*λa*(-1 + μt)) +δ*ηa*(θa**2)*(ϵ + θt - θt*μt)*τa -βa*θa*(ϵ + θt - θt*μt)*(δ*λa + γ*ηa*τa)) + βa*(ϵ + θt)*(βa*γ*θt*(-1 + μt) + δ*θa*(ϵ + θt - θt*μt))*ψt)) + θa*μa*(lt*βa*δ*ϵ*(ϵ + θt)*ρ -αt*(βa*δ*ϵ*λt*ρ - βa*γ*(ϵ + θt - θt*μt)*τa +δ*θa*(ϵ + θt - θt*μt)*τa) + βa*δ*ϵ*(ϵ + θt)*ρ*ψt)) +(lt*αa*βa*δ*ϵ*(ϵ + θt)*ρ - αa*αt*βa*δ*ϵ*λt*ρ + αa*αt*βa*γ*ϵ*τa -αa*αt*δ*ϵ*θa*τa + αa*αt*βa*γ*θt*τa - αa*αt*δ*θa*θt*τa -αa*αt*δ*ϵ*θa*μa*τa - αa*αt*δ*θa*θt*μa*τa - αa*αt*βa*γ*θt*μt*τa +αa*αt*δ*θa*θt*μt*τa + αa*αt*δ*θa*θt*μa*μt*τa + αt*βa*δ*ϵ*θa*μa*ρ*ψa +αa*βa*δ*(ϵ**2)*ρ*ψt + αa*βa*δ*ϵ*θt*ρ*ψt +e*(la*αt*βa*(βa*γ - δ*θa)*(ϵ + θt - θt*μt) -αt*βa*(βa*γ*θt*(-1 + μt) + δ*θa*(ϵ + θt - θt*μt))*ψa -αa*(ϵ + θt - θt*μt)*(αt*(βa*δ*λa + βa*γ*ηa*τa - 2*δ*ηa*θa*τa) + βa*δ*(ϵ + θt)*ψt))**2)))/(αa*αt*δ*(ϵ + θt - θt*μt)*(-(e*βa*γ) + βa*γ*ρ + δ*θa*(-1 + μa)*ρ))
  Ha <- -0.5*(e*la*αt*(βa**2)*γ*ϵ + e*la*αt*βa*δ*ϵ*θa + e*la*αt*(βa**2)*γ*θt + e*la*αt*βa*δ*θa*θt - e*αa*αt*βa*δ*ϵ*λa - e*αa*αt*βa*δ*θt*λa - e*la*αt*(βa**2)*γ*θt*μt - e*la*αt*βa*δ*θa*θt*μt + e*αa*αt*βa*δ*θt*λa*μt + lt*αa*βa*δ*(ϵ**2)*ρ + lt*αa*βa*δ*ϵ*θt*ρ - αa*αt*βa*δ*ϵ*λt*ρ + αa*αt*βa*γ*ϵ*τa - e*αa*αt*βa*γ*ϵ*ηa*τa + αa*αt*δ*ϵ*θa*τa + αa*αt*βa*γ*θt*τa - e*αa*αt*βa*γ*ηa*θt*τa + αa*αt*δ*θa*θt*τa - αa*αt*δ*ϵ*θa*μa*τa - αa*αt*δ*θa*θt*μa*τa - αa*αt*βa*γ*θt*μt*τa + e*αa*αt*βa*γ*ηa*θt*μt*τa - αa*αt*δ*θa*θt*μt*τa + αa*αt*δ*θa*θt*μa*μt*τa + e*αt*βa*δ*ϵ*θa*ψa + e*αt*(βa**2)*γ*θt*ψa + e*αt*βa*δ*θa*θt*ψa - e*αt*(βa**2)*γ*θt*μt*ψa - e*αt*βa*δ*θa*θt*μt*ψa - 2*αt*βa*δ*ϵ*θa*ρ*ψa + αt*βa*δ*ϵ*θa*μa*ρ*ψa - e*αa*βa*δ*(ϵ**2)*ψt - 2*e*αa*βa*δ*ϵ*θt*ψt - e*αa*βa*δ*(θt**2)*ψt + e*αa*βa*δ*ϵ*θt*μt*ψt + e*αa*βa*δ*(θt**2)*μt*ψt + αa*βa*δ*(ϵ**2)*ρ*ψt + αa*βa*δ*ϵ*θt*ρ*ψt - sqrt(4*αa*αt*δ*(αa*(ϵ + θt - θt*μt)*τa - βa*ϵ*ρ*ψa + e*(ϵ + θt - θt*μt)*(la*βa - αa*ηa*τa + βa*ψa))*(-(e*(lt*(βa**2)*γ*ϵ*(ϵ + θt) - αt*((βa**2)*γ*(ϵ*(λa + λt) - θt*λa*(-1 + μt)) +δ*ηa*(θa**2)*(ϵ + θt - θt*μt)*τa -βa*θa*(ϵ + θt - θt*μt)*(δ*λa + γ*ηa*τa)) + βa*(ϵ + θt)*(βa*γ*θt*(-1 + μt) + δ*θa*(ϵ + θt - θt*μt))*ψt)) + θa*μa*(lt*βa*δ*ϵ*(ϵ + θt)*ρ -αt*(βa*δ*ϵ*λt*ρ - βa*γ*(ϵ + θt - θt*μt)*τa +δ*θa*(ϵ + θt - θt*μt)*τa) + βa*δ*ϵ*(ϵ + θt)*ρ*ψt)) +(lt*αa*βa*δ*ϵ*(ϵ + θt)*ρ - αa*αt*βa*δ*ϵ*λt*ρ + αa*αt*βa*γ*ϵ*τa -αa*αt*δ*ϵ*θa*τa + αa*αt*βa*γ*θt*τa - αa*αt*δ*θa*θt*τa -αa*αt*δ*ϵ*θa*μa*τa - αa*αt*δ*θa*θt*μa*τa - αa*αt*βa*γ*θt*μt*τa +αa*αt*δ*θa*θt*μt*τa + αa*αt*δ*θa*θt*μa*μt*τa + αt*βa*δ*ϵ*θa*μa*ρ*ψa +αa*βa*δ*(ϵ**2)*ρ*ψt + αa*βa*δ*ϵ*θt*ρ*ψt +e*(la*αt*βa*(βa*γ - δ*θa)*(ϵ + θt - θt*μt) -αt*βa*(βa*γ*θt*(-1 + μt) + δ*θa*(ϵ + θt - θt*μt))*ψa -αa*(ϵ + θt - θt*μt)*(αt*(βa*δ*λa + βa*γ*ηa*τa - 2*δ*ηa*θa*τa) + βa*δ*(ϵ + θt)*ψt))**2)))/(αt*βa*δ*(αa*(ϵ + θt - θt*μt)*τa - βa*ϵ*ρ*ψa +e*(ϵ + θt - θt*μt)*(la*βa - αa*ηa*τa + βa*ψa)))
  Pa <- -0.5*(-(e*la*αt*(βa**2)*γ*ϵ*ρ) + e*la*αt*βa*δ*ϵ*θa*ρ - e*la*αt*(βa**2)*γ*θt*ρ + e*la*αt*βa*δ*θa*θt*ρ - e*αa*αt*βa*δ*ϵ*λa*ρ - e*αa*αt*βa*δ*θt*λa*ρ + e*la*αt*(βa**2)*γ*θt*μt*ρ - e*la*αt*βa*δ*θa*θt*μt*ρ + e*αa*αt*βa*δ*θt*λa*μt*ρ + lt*αa*βa*δ*(ϵ**2)*(ρ**2) + lt*αa*βa*δ*ϵ*θt*(ρ**2) - αa*αt*βa*δ*ϵ*λt*(ρ**2) - 2*e*αa*αt*βa*γ*ϵ*τa - 2*e*αa*αt*βa*γ*θt*τa + 2*e*αa*αt*βa*γ*θt*μt*τa + αa*αt*βa*γ*ϵ*ρ*τa + e*αa*αt*βa*γ*ϵ*ηa*ρ*τa - αa*αt*δ*ϵ*θa*ρ*τa + αa*αt*βa*γ*θt*ρ*τa + e*αa*αt*βa*γ*ηa*θt*ρ*τa - αa*αt*δ*θa*θt*ρ*τa + αa*αt*δ*ϵ*θa*μa*ρ*τa + αa*αt*δ*θa*θt*μa*ρ*τa - αa*αt*βa*γ*θt*μt*ρ*τa - e*αa*αt*βa*γ*ηa*θt*μt*ρ*τa + αa*αt*δ*θa*θt*μt*ρ*τa - αa*αt*δ*θa*θt*μa*μt*ρ*τa + e*αt*βa*δ*ϵ*θa*ρ*ψa - e*αt*(βa**2)*γ*θt*ρ*ψa + e*αt*βa*δ*θa*θt*ρ*ψa + e*αt*(βa**2)*γ*θt*μt*ρ*ψa - e*αt*βa*δ*θa*θt*μt*ρ*ψa - αt*βa*δ*ϵ*θa*μa*(ρ**2)*ψa - e*αa*βa*δ*(ϵ**2)*ρ*ψt - 2*e*αa*βa*δ*ϵ*θt*ρ*ψt - e*αa*βa*δ*(θt**2)*ρ*ψt + e*αa*βa*δ*ϵ*θt*μt*ρ*ψt + e*αa*βa*δ*(θt**2)*μt*ρ*ψt + αa*βa*δ*(ϵ**2)*(ρ**2)*ψt + αa*βa*δ*ϵ*θt*(ρ**2)*ψt + ρ*sqrt(4*αa*αt*δ*(αa*(ϵ + θt - θt*μt)*τa - βa*ϵ*ρ*ψa +e*(ϵ + θt - θt*μt)*(la*βa - αa*ηa*τa + βa*ψa))*(-(e*(lt*(βa**2)*γ*ϵ*(ϵ + θt) -αt*((βa**2)*γ*(ϵ*(λa + λt) - θt*λa*(-1 + μt)) +δ*ηa*(θa**2)*(ϵ + θt - θt*μt)*τa -βa*θa*(ϵ + θt - θt*μt)*(δ*λa + γ*ηa*τa)) +βa*(ϵ + θt)*(βa*γ*θt*(-1 + μt) + δ*θa*(ϵ + θt - θt*μt))*ψt)) +θa*μa*(lt*βa*δ*ϵ*(ϵ + θt)*ρ -αt*(βa*δ*ϵ*λt*ρ - βa*γ*(ϵ + θt - θt*μt)*τa +δ*θa*(ϵ + θt - θt*μt)*τa) + βa*δ*ϵ*(ϵ + θt)*ρ*ψt)) +(lt*αa*βa*δ*ϵ*(ϵ + θt)*ρ - αa*αt*βa*δ*ϵ*λt*ρ + αa*αt*βa*γ*ϵ*τa - αa*αt*δ*ϵ*θa*τa + αa*αt*βa*γ*θt*τa - αa*αt*δ*θa*θt*τa - αa*αt*δ*ϵ*θa*μa*τa - αa*αt*δ*θa*θt*μa*τa - αa*αt*βa*γ*θt*μt*τa + αa*αt*δ*θa*θt*μt*τa + αa*αt*δ*θa*θt*μa*μt*τa + αt*βa*δ*ϵ*θa*μa*ρ*ψa + αa*βa*δ*(ϵ**2)*ρ*ψt + αa*βa*δ*ϵ*θt*ρ*ψt + e*(la*αt*βa*(βa*γ - δ*θa)*(ϵ + θt - θt*μt) -αt*βa*(βa*γ*θt*(-1 + μt) + δ*θa*(ϵ + θt - θt*μt))*ψa -αa*(ϵ + θt - θt*μt)*(αt*(βa*δ*λa + βa*γ*ηa*τa - 2*δ*ηa*θa*τa) + βa*δ*(ϵ + θt)*ψt))**2)))/(e*αa*αt*βa*(ϵ + θt - θt*μt)*(e*βa*γ - (βa*γ + δ*θa*(-1 + μa))*ρ))
  Na <- -0.5*(-(e*la*αt*(βa**2)*γ*ϵ) + e*la*αt*βa*δ*ϵ*θa - e*la*αt*(βa**2)*γ*θt + e*la*αt*βa*δ*θa*θt + e*αa*αt*βa*δ*ϵ*λa + e*αa*αt*βa*δ*θt*λa + e*la*αt*(βa**2)*γ*θt*μt - e*la*αt*βa*δ*θa*θt*μt - e*αa*αt*βa*δ*θt*λa*μt - lt*αa*βa*δ*(ϵ**2)*ρ - lt*αa*βa*δ*ϵ*θt*ρ + αa*αt*βa*δ*ϵ*λt*ρ - αa*αt*βa*γ*ϵ*τa + e*αa*αt*βa*γ*ϵ*ηa*τa + αa*αt*δ*ϵ*θa*τa - 2*e*αa*αt*δ*ϵ*ηa*θa*τa - αa*αt*βa*γ*θt*τa + e*αa*αt*βa*γ*ηa*θt*τa + αa*αt*δ*θa*θt*τa - 2*e*αa*αt*δ*ηa*θa*θt*τa + αa*αt*δ*ϵ*θa*μa*τa + αa*αt*δ*θa*θt*μa*τa + αa*αt*βa*γ*θt*μt*τa - e*αa*αt*βa*γ*ηa*θt*μt*τa - αa*αt*δ*θa*θt*μt*τa + 2*e*αa*αt*δ*ηa*θa*θt*μt*τa - αa*αt*δ*θa*θt*μa*μt*τa + e*αt*βa*δ*ϵ*θa*ψa - e*αt*(βa**2)*γ*θt*ψa + e*αt*βa*δ*θa*θt*ψa + e*αt*(βa**2)*γ*θt*μt*ψa - e*αt*βa*δ*θa*θt*μt*ψa - αt*βa*δ*ϵ*θa*μa*ρ*ψa + e*αa*βa*δ*(ϵ**2)*ψt + 2*e*αa*βa*δ*ϵ*θt*ψt + e*αa*βa*δ*(θt**2)*ψt - e*αa*βa*δ*ϵ*θt*μt*ψt - e*αa*βa*δ*(θt**2)*μt*ψt - αa*βa*δ*(ϵ**2)*ρ*ψt - αa*βa*δ*ϵ*θt*ρ*ψt + sqrt(4*αa*αt*δ*(αa*(ϵ + θt - θt*μt)*τa - βa*ϵ*ρ*ψa + e*(ϵ + θt - θt*μt)*(la*βa - αa*ηa*τa + βa*ψa))*(-(e*(lt*(βa**2)*γ*ϵ*(ϵ + θt) - αt*((βa**2)*γ*(ϵ*(λa + λt) - θt*λa*(-1 + μt)) +δ*ηa*(θa**2)*(ϵ + θt - θt*μt)*τa -βa*θa*(ϵ + θt - θt*μt)*(δ*λa + γ*ηa*τa)) + βa*(ϵ + θt)*(βa*γ*θt*(-1 + μt) + δ*θa*(ϵ + θt - θt*μt))*ψt)) + θa*μa*(lt*βa*δ*ϵ*(ϵ + θt)*ρ -αt*(βa*δ*ϵ*λt*ρ - βa*γ*(ϵ + θt - θt*μt)*τa +δ*θa*(ϵ + θt - θt*μt)*τa) + βa*δ*ϵ*(ϵ + θt)*ρ*ψt)) +(lt*αa*βa*δ*ϵ*(ϵ + θt)*ρ - αa*αt*βa*δ*ϵ*λt*ρ + αa*αt*βa*γ*ϵ*τa -αa*αt*δ*ϵ*θa*τa + αa*αt*βa*γ*θt*τa - αa*αt*δ*θa*θt*τa -αa*αt*δ*ϵ*θa*μa*τa - αa*αt*δ*θa*θt*μa*τa - αa*αt*βa*γ*θt*μt*τa +αa*αt*δ*θa*θt*μt*τa + αa*αt*δ*θa*θt*μa*μt*τa + αt*βa*δ*ϵ*θa*μa*ρ*ψa +αa*βa*δ*(ϵ**2)*ρ*ψt + αa*βa*δ*ϵ*θt*ρ*ψt +e*(la*αt*βa*(βa*γ - δ*θa)*(ϵ + θt - θt*μt) -αt*βa*(βa*γ*θt*(-1 + μt) + δ*θa*(ϵ + θt - θt*μt))*ψa -αa*(ϵ + θt - θt*μt)*(αt*(βa*δ*λa + βa*γ*ηa*τa - 2*δ*ηa*θa*τa) + βa*δ*(ϵ + θt)*ψt))**2)))/(αa*αt*δ*(-(αa*(ϵ + θt - θt*μt)*τa) + βa*ϵ*ρ*ψa -e*(ϵ + θt - θt*μt)*(la*βa - αa*ηa*τa + βa*ψa)))
  ha_prod <- βa*e*Ha*Pa+δ*ρ*Lt*Ha
  pt_prod <- αt*Pt*Nt
  pa_prod <- αa*Pa*Na
  pt_flux <- ϵ*Pt
  nt_flux <- ψt*Nt
  na_flux <- ψa*Na
  
  new_row <- c(αt,
                αa, 
                βa, 
                e,
                ρ,
                δ, 
                γ,
                ϵ,
                ηa, 
                θt, 
                θa,
                λt, 
                λa, 
                μt,
                μa,
                τa, 
                ψt,
                ψa,
                lt,
                la,
                Pt,
                Nt,
                Lt,
                Ha,
                Pa,
                Na,
                ha_prod,
                pt_prod,
                pa_prod,
                pt_flux,
                nt_flux,
                na_flux)
  
  param_simulations[i,] <- new_row
  
}
toc()

# filter for feasibility
tic("filter for feasibility")
param_simulations_feasible <- param_simulations %>% 
  filter(Pt > 0 & 
           Nt > 0 & 
           Lt > 0 & 
           Ha > 0 & 
           Pa > 0 & 
           Na > 0)
toc()

# calculate stability
tic("filter for stability")
stability <- NULL

for (i in 1:nrow(param_simulations)) {
  # fetch parameter values from the equilibrium simulations df
  αt = param_simulations_feasible$αt[i]
  αa = param_simulations_feasible$αa[i]
  βa = param_simulations_feasible$βa[i]
  e = param_simulations_feasible$e[i]
  ρ = param_simulations_feasible$ρ[i]
  δ = param_simulations_feasible$δ[i]
  γ = param_simulations_feasible$γ[i]
  ϵ = param_simulations_feasible$ϵ[i]
  ηa = param_simulations_feasible$ηa[i]
  θt = param_simulations_feasible$θt[i]
  θa = param_simulations_feasible$θa[i]
  λt = param_simulations_feasible$λt[i]
  λa = param_simulations_feasible$λa[i]
  μt = param_simulations_feasible$μt[i]
  μa = param_simulations_feasible$μa[i]
  τa = param_simulations_feasible$τa[i]
  ψt = param_simulations_feasible$ψt[i]
  ψa = param_simulations_feasible$ψa[i]
  lt = param_simulations_feasible$lt[i]
  la = param_simulations_feasible$la[i]
  
  # Put in the solutions to the equilibrium values here:
  dPt = param_simulations_feasible$Pt[i]
  dNt = param_simulations_feasible$Nt[i]
  dLt = param_simulations_feasible$Lt[i]
  dHa = param_simulations_feasible$Ha[i]
  dPa = param_simulations_feasible$Pa[i]
  dNa = param_simulations_feasible$Na[i]
  
  # create the Jacobian from Mathematica
  jacob <- rbind(c(dNt*αt-ϵ-θt, dPt*αt, 0, 0, 0, 0), 
                 c(-(dNt*αt)+θt*μt, -lt-dPt*αt-ψt, 0, 0, 0, ψa), 
                 c(ϵ, 0, -γ-dHa*δ, -(dLt*δ), 0, 0), 
                 c(0, 0, dHa*δ*ρ, e*dPa*βa+dLt*δ*ρ-τa, e*dHa*βa, 0), 
                 c(0, 0, 0, -(dPa*βa), dNa*αa-dHa*βa-θa, dPa*αa), 
                 c(0, ψt, γ, ηa*τa, -(dNa*αa)+θa*μa, -la-dPa*αa-ψa))
  
  # Check to see if the real part of the leading eigenvalue is less than 0 or not - if it is than the system is stable.
  if (max(Re(base::eigen(jacob)$values)) < 0 ){ 
    stable <- "stable"
  } else {
    stable <- "unstable"
  }
  
  stability <- rbind(stability, data.frame(αt,
                                            αa, 
                                            βa, 
                                            e,
                                            ρ,
                                            δ, 
                                            γ,
                                            ϵ,
                                            ηa, 
                                            θt, 
                                            θa,
                                            λt, 
                                            λa, 
                                            μt,
                                            μa, 
                                            τa, 
                                            ψt,
                                            ψa,
                                            lt,
                                            la,
                                            dPt, 
                                            dNt,
                                            dLt,
                                            dHa, 
                                            dPa, 
                                            dNa, 
                                            EiV1 = Re(eigen(jacob)$values[1]),
                                            EiV2 = Re(eigen(jacob)$values[2]),
                                            EiV3 = Re(eigen(jacob)$values[3]),
                                            EiV4 = Re(eigen(jacob)$values[4]),
                                            EiV5 = Re(eigen(jacob)$values[5]),
                                            EiV6 = Re(eigen(jacob)$values[6]),
                                            maxEv = max(Re(base::eigen(jacob)$values)),
                                            stable = stable))
  
}

# join stability and param simulations together and filter out rows that produce unstable or negative equilibrium values
param_simulations_stability_test <- param_simulations_feasible %>% 
  left_join(stability %>% select(-dPt, 
                                 -dNt,
                                 -dLt,
                                 -dHa, 
                                 -dPa, 
                                 -dNa, 
                                 -EiV1, 
                                 -EiV2, 
                                 -EiV3, 
                                 -EiV4, 
                                 -EiV5,
                                 -EiV6,
                                 -maxEv), by = c("αt",
                                                 "αa",
                                                 "βa",
                                                 "e",
                                                 "ρ",
                                                 "δ",
                                                 "γ",
                                                 "ϵ",
                                                 "ηa",
                                                 "θt",
                                                 "θa",
                                                 "λt",
                                                 "λa",
                                                 "μt",
                                                 "μa",
                                                 "τa",
                                                 "ψt",
                                                 "ψa",
                                                 "lt",
                                                 "la"))

param_simulations_stable <- param_simulations_stability_test %>% 
  filter(stable == "stable")
toc()

# export file with all feasible and stable parameter combinations
write.csv(param_simulations_stable, "output/param_simulations_stable_10000.csv", row.names=FALSE)

# first make a function that takes in one row and generates a dataset for the sensitivity test
expand_row <- function(master_df,
                       row_number,
                       param_number,
                       min_param_val,
                       max_param_val,
                       increase_val,
                       expand_length){
  # master_df: dataframe with original parameters to be expanded
  # row_number: row number in master_df to expand
  # param_number: column number of the parameter to be changed for sensitivity test
  expanded_data<-data.frame(αt=master_df$αt[row_number], 
                             αa=master_df$αa[row_number],
                             βa=master_df$βa[row_number],
                             e=master_df$e[row_number],
                             δ=master_df$δ[row_number],
                             γ=master_df$γ[row_number],
                             ϵ=master_df$ϵ[row_number],
                             ηa=master_df$ηa[row_number],
                             θt=master_df$θt[row_number],
                             θa=master_df$θa[row_number],
                             λt=master_df$λt[row_number],
                             λa=master_df$λa[row_number],
                             μt=master_df$μt[row_number],
                             μa=master_df$μa[row_number], 
                             ρ=master_df$ρ[row_number],
                             τa=master_df$τa[row_number], 
                             ψt=master_df$ψt[row_number],
                             ψa=master_df$ψa[row_number],
                             lt=master_df$lt[row_number],
                             la=master_df$la[row_number]) %>% 
    slice(rep(1:n(), each =  expand_length))
  
  # create a vector of the new parameter values
  param_vector <- seq(min_param_val,max_param_val,increase_val)
  
  # now change the desired column to this vector
  expanded_data[,param_number] = param_vector
  
  return(expanded_data)
  
}

# now make a larger function to iterate through simulated parameters, make a 
# list of expanded data frames, and append together

simulate_disturbance <- function(master_df,
                             param_number, 
                             min_param_val,
                             max_param_val, 
                             increase_val,
                             expand_length,
                             param_number2, 
                             min_param_val2,
                             max_param_val2, 
                             increase_val2,
                             expand_length2){
  # Param number: number of the parameter in params (e.g., alpha = 1, gamma=2)
  # max_param_val: maximum value for the parameter to be tested to
  # increase_val: value that the parameter will be increased by for each iteration
  
  # make an empty list to add the expanded data frames onto
  expanded_data_main <- NULL
  
  for (i in 1:nrow(master_df)){ 
    
    expanded_df <- expand_row(master_df,
                              row_number=i,
                              param_number,
                              min_param_val,
                              max_param_val,
                              increase_val,
                              expand_length)
    
    expanded_data_main <- expanded_data_main %>% 
      bind_rows(expanded_df)
  }
  
  
  # repeat for a second variable
  expanded_data_main2 <- NULL
  for (i in 1:nrow(expanded_data_main)){ 
    
    expanded_df2 <- expand_row(expanded_data_main,
                               row_number=i,
                               param_number2,
                               min_param_val2,
                               max_param_val2,
                               increase_val2,
                               expand_length2)
    
    expanded_data_main2 <- expanded_data_main2 %>% 
      bind_rows(expanded_df2)
  }
  
  
  
  # recalculate equilibrium value
  sensitivity_params <- expanded_data_main2 %>% mutate(Pt = -0.5*(-2*e*la*lt*αa*βa*δ*(ϵ**2) - 4*e*la*lt*αa*βa*δ*ϵ*θt -  2*e*la*lt*αa*βa*δ*(θt**2) + 2*e*la*αa*αt*βa*δ*ϵ*λt +  
                                                                    2*e*la*αa*αt*βa*δ*θt*λt + 2*e*la*lt*αa*βa*δ*ϵ*θt*μt +  2*e*la*lt*αa*βa*δ*(θt**2)*μt - 2*e*la*αa*αt*βa*δ*θt*λt*μt -  
                                                                    2*lt*(αa**2)*δ*(ϵ**2)*τa + 2*e*lt*(αa**2)*δ*(ϵ**2)*ηa*τa -  4*lt*(αa**2)*δ*ϵ*θt*τa + 4*e*lt*(αa**2)*δ*ϵ*ηa*θt*τa -  
                                                                    2*lt*(αa**2)*δ*(θt**2)*τa + 2*e*lt*(αa**2)*δ*ηa*(θt**2)*τa +  2*(αa**2)*αt*δ*ϵ*λt*τa - 2*e*(αa**2)*αt*δ*ϵ*ηa*λt*τa +  
                                                                    2*(αa**2)*αt*δ*θt*λt*τa - 2*e*(αa**2)*αt*δ*ηa*θt*λt*τa +  2*lt*(αa**2)*δ*ϵ*θt*μt*τa - 2*e*lt*(αa**2)*δ*ϵ*ηa*θt*μt*τa +  
                                                                    2*lt*(αa**2)*δ*(θt**2)*μt*τa -  2*e*lt*(αa**2)*δ*ηa*(θt**2)*μt*τa - 2*(αa**2)*αt*δ*θt*λt*μt*τa +  2*e*(αa**2)*αt*δ*ηa*θt*λt*μt*τa - 
                                                                    e*la*αt*(βa**2)*γ*ϵ*ψa -  2*e*lt*αa*βa*δ*(ϵ**2)*ψa + e*la*αt*βa*δ*ϵ*θa*ψa -  e*la*αt*(βa**2)*γ*θt*ψa - 4*e*lt*αa*βa*δ*ϵ*θt*ψa +  
                                                                    e*la*αt*βa*δ*θa*θt*ψa - 2*e*lt*αa*βa*δ*(θt**2)*ψa +  e*αa*αt*βa*δ*ϵ*λa*ψa + e*αa*αt*βa*δ*θt*λa*ψa + 2*e*αa*αt*βa*δ*ϵ*λt*ψa +  
                                                                    2*e*αa*αt*βa*δ*θt*λt*ψa + e*la*αt*(βa**2)*γ*θt*μt*ψa +  2*e*lt*αa*βa*δ*ϵ*θt*μt*ψa - e*la*αt*βa*δ*θa*θt*μt*ψa +  
                                                                    2*e*lt*αa*βa*δ*(θt**2)*μt*ψa - e*αa*αt*βa*δ*θt*λa*μt*ψa -  2*e*αa*αt*βa*δ*θt*λt*μt*ψa + lt*αa*βa*δ*(ϵ**2)*ρ*ψa +  
                                                                    lt*αa*βa*δ*ϵ*θt*ρ*ψa - αa*αt*βa*δ*ϵ*λt*ρ*ψa - αa*αt*βa*γ*ϵ*τa*ψa +  e*αa*αt*βa*γ*ϵ*ηa*τa*ψa + αa*αt*δ*ϵ*θa*τa*ψa - 
                                                                    2*e*αa*αt*δ*ϵ*ηa*θa*τa*ψa -  αa*αt*βa*γ*θt*τa*ψa + e*αa*αt*βa*γ*ηa*θt*τa*ψa + αa*αt*δ*θa*θt*τa*ψa -  
                                                                    2*e*αa*αt*δ*ηa*θa*θt*τa*ψa + αa*αt*δ*ϵ*θa*μa*τa*ψa + αa*αt*δ*θa*θt*μa*τa*ψa +  αa*αt*βa*γ*θt*μt*τa*ψa - 
                                                                    e*αa*αt*βa*γ*ηa*θt*μt*τa*ψa -  αa*αt*δ*θa*θt*μt*τa*ψa + 2*e*αa*αt*δ*ηa*θa*θt*μt*τa*ψa -  αa*αt*δ*θa*θt*μa*μt*τa*ψa + 
                                                                    e*αt*βa*δ*ϵ*θa*(ψa**2) -  e*αt*(βa**2)*γ*θt*(ψa**2) + e*αt*βa*δ*θa*θt*(ψa**2) +  e*αt*(βa**2)*γ*θt*μt*(ψa**2) - 
                                                                    e*αt*βa*δ*θa*θt*μt*(ψa**2) -  αt*βa*δ*ϵ*θa*μa*ρ*(ψa**2) - 2*e*la*αa*βa*δ*(ϵ**2)*ψt -  4*e*la*αa*βa*δ*ϵ*θt*ψt - 
                                                                    2*e*la*αa*βa*δ*(θt**2)*ψt +  2*e*la*αa*βa*δ*ϵ*θt*μt*ψt + 2*e*la*αa*βa*δ*(θt**2)*μt*ψt -  2*(αa**2)*δ*(ϵ**2)*τa*ψt + 
                                                                    2*e*(αa**2)*δ*(ϵ**2)*ηa*τa*ψt -  4*(αa**2)*δ*ϵ*θt*τa*ψt + 4*e*(αa**2)*δ*ϵ*ηa*θt*τa*ψt -  2*(αa**2)*δ*(θt**2)*τa*ψt + 
                                                                    2*e*(αa**2)*δ*ηa*(θt**2)*τa*ψt +  2*(αa**2)*δ*ϵ*θt*μt*τa*ψt - 2*e*(αa**2)*δ*ϵ*ηa*θt*μt*τa*ψt +  2*(αa**2)*δ*(θt**2)*μt*τa*ψt -  
                                                                    2*e*(αa**2)*δ*ηa*(θt**2)*μt*τa*ψt - e*αa*βa*δ*(ϵ**2)*ψa*ψt -  2*e*αa*βa*δ*ϵ*θt*ψa*ψt - e*αa*βa*δ*(θt**2)*ψa*ψt +  
                                                                    e*αa*βa*δ*ϵ*θt*μt*ψa*ψt + e*αa*βa*δ*(θt**2)*μt*ψa*ψt +  αa*βa*δ*(ϵ**2)*ρ*ψa*ψt + αa*βa*δ*ϵ*θt*ρ*ψa*ψt +  
                                                                    ψa*sqrt(4*αa*αt*δ*(αa*(ϵ + θt - θt*μt)*τa - βa*ϵ*ρ*ψa + e*(ϵ + θt - θt*μt)*(la*βa - αa*ηa*τa + βa*ψa))* 
                                                                              (-(e*(lt*(βa**2)*γ*ϵ*(ϵ + θt) -  αt*((βa**2)*γ*(ϵ*(λa + λt) - θt*λa*(-1 + μt)) +  δ*ηa*(θa**2)*
                                                                                                                     (ϵ + θt - θt*μt)*τa -  βa*θa*(ϵ + θt - θt*μt)*(δ*λa + γ*ηa*τa)) + 
                                                                                      βa*(ϵ + θt)*(βa*γ*θt*(-1 + μt) + δ*θa*(ϵ + θt - θt*μt))*ψt)) + 
                                                                                 θa*μa*(lt*βa*δ*ϵ*(ϵ + θt)*ρ -αt*(βa*δ*ϵ*λt*ρ - βa*γ*(ϵ + θt - θt*μt)*τa +δ*θa*(ϵ + θt - θt*μt)*τa) + 
                                                                                          βa*δ*ϵ*(ϵ + θt)*ρ*ψt)) + (lt*αa*βa*δ*ϵ*(ϵ + θt)*ρ - αa*αt*βa*δ*ϵ*λt*ρ + αa*αt*βa*γ*ϵ*τa -αa*αt*δ*ϵ*θa*τa + 
                                                                                                                      αa*αt*βa*γ*θt*τa - αa*αt*δ*θa*θt*τa -αa*αt*δ*ϵ*θa*μa*τa - αa*αt*δ*θa*θt*μa*τa - 
                                                                                                                      αa*αt*βa*γ*θt*μt*τa +αa*αt*δ*θa*θt*μt*τa + αa*αt*δ*θa*θt*μa*μt*τa + αt*βa*δ*ϵ*θa*μa*ρ*ψa +
                                                                                                                      αa*βa*δ*(ϵ**2)*ρ*ψt + αa*βa*δ*ϵ*θt*ρ*ψt +
                                                                                                                      e*(la*αt*βa*(βa*γ - δ*θa)*(ϵ + θt - θt*μt) -   αt*βa*(βa*γ*θt*(-1 + μt) + δ*θa*(ϵ + θt - θt*μt))*ψa -   
                                                                                                                           αa*(ϵ + θt - θt*μt)*   (αt*(βa*δ*λa + βa*γ*ηa*τa - 2*δ*ηa*θa*τa) + βa*δ*(ϵ + θt)*ψt))**2)))/ 
                                                         (αa*αt*δ*(ϵ + θt - θt*μt)*(-(αa*(ϵ + θt - θt*μt)*τa) + βa*ϵ*ρ*ψa -e*(ϵ + θt - θt*μt)*(la*βa - αa*ηa*τa + βa*ψa))),
                                                       Nt = (ϵ + θt)/αt,
                                                       Lt = -0.5*(-(e*la*αt*(βa**2)*γ*ϵ) + e*la*αt*βa*δ*ϵ*θa - e*la*αt*(βa**2)*γ*θt + e*la*αt*βa*δ*θa*θt - e*αa*αt*βa*δ*ϵ*λa - e*αa*αt*βa*δ*θt*λa + e*la*αt*(βa**2)*γ*θt*μt - e*la*αt*βa*δ*θa*θt*μt + e*αa*αt*βa*δ*θt*λa*μt + lt*αa*βa*δ*(ϵ**2)*ρ + lt*αa*βa*δ*ϵ*θt*ρ - αa*αt*βa*δ*ϵ*λt*ρ - αa*αt*βa*γ*ϵ*τa + e*αa*αt*βa*γ*ϵ*ηa*τa + αa*αt*δ*ϵ*θa*τa - αa*αt*βa*γ*θt*τa + e*αa*αt*βa*γ*ηa*θt*τa + αa*αt*δ*θa*θt*τa - αa*αt*δ*ϵ*θa*μa*τa - αa*αt*δ*θa*θt*μa*τa + αa*αt*βa*γ*θt*μt*τa - e*αa*αt*βa*γ*ηa*θt*μt*τa - αa*αt*δ*θa*θt*μt*τa + αa*αt*δ*θa*θt*μa*μt*τa + e*αt*βa*δ*ϵ*θa*ψa - e*αt*(βa**2)*γ*θt*ψa + e*αt*βa*δ*θa*θt*ψa + e*αt*(βa**2)*γ*θt*μt*ψa - e*αt*βa*δ*θa*θt*μt*ψa - αt*βa*δ*ϵ*θa*μa*ρ*ψa - e*αa*βa*δ*(ϵ**2)*ψt - 2*e*αa*βa*δ*ϵ*θt*ψt - e*αa*βa*δ*(θt**2)*ψt + e*αa*βa*δ*ϵ*θt*μt*ψt + e*αa*βa*δ*(θt**2)*μt*ψt + αa*βa*δ*(ϵ**2)*ρ*ψt + αa*βa*δ*ϵ*θt*ρ*ψt + sqrt(4*αa*αt*δ*(αa*(ϵ + θt - θt*μt)*τa - βa*ϵ*ρ*ψa + e*(ϵ + θt - θt*μt)*(la*βa - αa*ηa*τa + βa*ψa))*(-(e*(lt*(βa**2)*γ*ϵ*(ϵ + θt) - αt*((βa**2)*γ*(ϵ*(λa + λt) - θt*λa*(-1 + μt)) +δ*ηa*(θa**2)*(ϵ + θt - θt*μt)*τa -βa*θa*(ϵ + θt - θt*μt)*(δ*λa + γ*ηa*τa)) + βa*(ϵ + θt)*(βa*γ*θt*(-1 + μt) + δ*θa*(ϵ + θt - θt*μt))*ψt)) + θa*μa*(lt*βa*δ*ϵ*(ϵ + θt)*ρ -αt*(βa*δ*ϵ*λt*ρ - βa*γ*(ϵ + θt - θt*μt)*τa +δ*θa*(ϵ + θt - θt*μt)*τa) + βa*δ*ϵ*(ϵ + θt)*ρ*ψt)) +(lt*αa*βa*δ*ϵ*(ϵ + θt)*ρ - αa*αt*βa*δ*ϵ*λt*ρ + αa*αt*βa*γ*ϵ*τa -αa*αt*δ*ϵ*θa*τa + αa*αt*βa*γ*θt*τa - αa*αt*δ*θa*θt*τa -αa*αt*δ*ϵ*θa*μa*τa - αa*αt*δ*θa*θt*μa*τa - αa*αt*βa*γ*θt*μt*τa +αa*αt*δ*θa*θt*μt*τa + αa*αt*δ*θa*θt*μa*μt*τa + αt*βa*δ*ϵ*θa*μa*ρ*ψa +αa*βa*δ*(ϵ**2)*ρ*ψt + αa*βa*δ*ϵ*θt*ρ*ψt +e*(la*αt*βa*(βa*γ - δ*θa)*(ϵ + θt - θt*μt) -αt*βa*(βa*γ*θt*(-1 + μt) + δ*θa*(ϵ + θt - θt*μt))*ψa -αa*(ϵ + θt - θt*μt)*(αt*(βa*δ*λa + βa*γ*ηa*τa - 2*δ*ηa*θa*τa) + βa*δ*(ϵ + θt)*ψt))**2)))/(αa*αt*δ*(ϵ + θt - θt*μt)*(-(e*βa*γ) + βa*γ*ρ + δ*θa*(-1 + μa)*ρ)),
                                                       Ha = -0.5*(e*la*αt*(βa**2)*γ*ϵ + e*la*αt*βa*δ*ϵ*θa + e*la*αt*(βa**2)*γ*θt + e*la*αt*βa*δ*θa*θt - e*αa*αt*βa*δ*ϵ*λa - e*αa*αt*βa*δ*θt*λa - e*la*αt*(βa**2)*γ*θt*μt - e*la*αt*βa*δ*θa*θt*μt + e*αa*αt*βa*δ*θt*λa*μt + lt*αa*βa*δ*(ϵ**2)*ρ + lt*αa*βa*δ*ϵ*θt*ρ - αa*αt*βa*δ*ϵ*λt*ρ + αa*αt*βa*γ*ϵ*τa - e*αa*αt*βa*γ*ϵ*ηa*τa + αa*αt*δ*ϵ*θa*τa + αa*αt*βa*γ*θt*τa - e*αa*αt*βa*γ*ηa*θt*τa + αa*αt*δ*θa*θt*τa - αa*αt*δ*ϵ*θa*μa*τa - αa*αt*δ*θa*θt*μa*τa - αa*αt*βa*γ*θt*μt*τa + e*αa*αt*βa*γ*ηa*θt*μt*τa - αa*αt*δ*θa*θt*μt*τa + αa*αt*δ*θa*θt*μa*μt*τa + e*αt*βa*δ*ϵ*θa*ψa + e*αt*(βa**2)*γ*θt*ψa + e*αt*βa*δ*θa*θt*ψa - e*αt*(βa**2)*γ*θt*μt*ψa - e*αt*βa*δ*θa*θt*μt*ψa - 2*αt*βa*δ*ϵ*θa*ρ*ψa + αt*βa*δ*ϵ*θa*μa*ρ*ψa - e*αa*βa*δ*(ϵ**2)*ψt - 2*e*αa*βa*δ*ϵ*θt*ψt - e*αa*βa*δ*(θt**2)*ψt + e*αa*βa*δ*ϵ*θt*μt*ψt + e*αa*βa*δ*(θt**2)*μt*ψt + αa*βa*δ*(ϵ**2)*ρ*ψt + αa*βa*δ*ϵ*θt*ρ*ψt - sqrt(4*αa*αt*δ*(αa*(ϵ + θt - θt*μt)*τa - βa*ϵ*ρ*ψa + e*(ϵ + θt - θt*μt)*(la*βa - αa*ηa*τa + βa*ψa))*(-(e*(lt*(βa**2)*γ*ϵ*(ϵ + θt) - αt*((βa**2)*γ*(ϵ*(λa + λt) - θt*λa*(-1 + μt)) +δ*ηa*(θa**2)*(ϵ + θt - θt*μt)*τa -βa*θa*(ϵ + θt - θt*μt)*(δ*λa + γ*ηa*τa)) + βa*(ϵ + θt)*(βa*γ*θt*(-1 + μt) + δ*θa*(ϵ + θt - θt*μt))*ψt)) + θa*μa*(lt*βa*δ*ϵ*(ϵ + θt)*ρ -αt*(βa*δ*ϵ*λt*ρ - βa*γ*(ϵ + θt - θt*μt)*τa +δ*θa*(ϵ + θt - θt*μt)*τa) + βa*δ*ϵ*(ϵ + θt)*ρ*ψt)) +(lt*αa*βa*δ*ϵ*(ϵ + θt)*ρ - αa*αt*βa*δ*ϵ*λt*ρ + αa*αt*βa*γ*ϵ*τa -αa*αt*δ*ϵ*θa*τa + αa*αt*βa*γ*θt*τa - αa*αt*δ*θa*θt*τa -αa*αt*δ*ϵ*θa*μa*τa - αa*αt*δ*θa*θt*μa*τa - αa*αt*βa*γ*θt*μt*τa +αa*αt*δ*θa*θt*μt*τa + αa*αt*δ*θa*θt*μa*μt*τa + αt*βa*δ*ϵ*θa*μa*ρ*ψa +αa*βa*δ*(ϵ**2)*ρ*ψt + αa*βa*δ*ϵ*θt*ρ*ψt +e*(la*αt*βa*(βa*γ - δ*θa)*(ϵ + θt - θt*μt) -αt*βa*(βa*γ*θt*(-1 + μt) + δ*θa*(ϵ + θt - θt*μt))*ψa -αa*(ϵ + θt - θt*μt)*(αt*(βa*δ*λa + βa*γ*ηa*τa - 2*δ*ηa*θa*τa) + βa*δ*(ϵ + θt)*ψt))**2)))/(αt*βa*δ*(αa*(ϵ + θt - θt*μt)*τa - βa*ϵ*ρ*ψa +e*(ϵ + θt - θt*μt)*(la*βa - αa*ηa*τa + βa*ψa))),
                                                       Pa = -0.5*(-(e*la*αt*(βa**2)*γ*ϵ*ρ) + e*la*αt*βa*δ*ϵ*θa*ρ - e*la*αt*(βa**2)*γ*θt*ρ + e*la*αt*βa*δ*θa*θt*ρ - e*αa*αt*βa*δ*ϵ*λa*ρ - e*αa*αt*βa*δ*θt*λa*ρ + e*la*αt*(βa**2)*γ*θt*μt*ρ - e*la*αt*βa*δ*θa*θt*μt*ρ + e*αa*αt*βa*δ*θt*λa*μt*ρ + lt*αa*βa*δ*(ϵ**2)*(ρ**2) + lt*αa*βa*δ*ϵ*θt*(ρ**2) - αa*αt*βa*δ*ϵ*λt*(ρ**2) - 2*e*αa*αt*βa*γ*ϵ*τa - 2*e*αa*αt*βa*γ*θt*τa + 2*e*αa*αt*βa*γ*θt*μt*τa + αa*αt*βa*γ*ϵ*ρ*τa + e*αa*αt*βa*γ*ϵ*ηa*ρ*τa - αa*αt*δ*ϵ*θa*ρ*τa + αa*αt*βa*γ*θt*ρ*τa + e*αa*αt*βa*γ*ηa*θt*ρ*τa - αa*αt*δ*θa*θt*ρ*τa + αa*αt*δ*ϵ*θa*μa*ρ*τa + αa*αt*δ*θa*θt*μa*ρ*τa - αa*αt*βa*γ*θt*μt*ρ*τa - e*αa*αt*βa*γ*ηa*θt*μt*ρ*τa + αa*αt*δ*θa*θt*μt*ρ*τa - αa*αt*δ*θa*θt*μa*μt*ρ*τa + e*αt*βa*δ*ϵ*θa*ρ*ψa - e*αt*(βa**2)*γ*θt*ρ*ψa + e*αt*βa*δ*θa*θt*ρ*ψa + e*αt*(βa**2)*γ*θt*μt*ρ*ψa - e*αt*βa*δ*θa*θt*μt*ρ*ψa - αt*βa*δ*ϵ*θa*μa*(ρ**2)*ψa - e*αa*βa*δ*(ϵ**2)*ρ*ψt - 2*e*αa*βa*δ*ϵ*θt*ρ*ψt - e*αa*βa*δ*(θt**2)*ρ*ψt + e*αa*βa*δ*ϵ*θt*μt*ρ*ψt + e*αa*βa*δ*(θt**2)*μt*ρ*ψt + αa*βa*δ*(ϵ**2)*(ρ**2)*ψt + αa*βa*δ*ϵ*θt*(ρ**2)*ψt + ρ*sqrt(4*αa*αt*δ*(αa*(ϵ + θt - θt*μt)*τa - βa*ϵ*ρ*ψa +e*(ϵ + θt - θt*μt)*(la*βa - αa*ηa*τa + βa*ψa))*(-(e*(lt*(βa**2)*γ*ϵ*(ϵ + θt) -αt*((βa**2)*γ*(ϵ*(λa + λt) - θt*λa*(-1 + μt)) +δ*ηa*(θa**2)*(ϵ + θt - θt*μt)*τa -βa*θa*(ϵ + θt - θt*μt)*(δ*λa + γ*ηa*τa)) +βa*(ϵ + θt)*(βa*γ*θt*(-1 + μt) + δ*θa*(ϵ + θt - θt*μt))*ψt)) +θa*μa*(lt*βa*δ*ϵ*(ϵ + θt)*ρ -αt*(βa*δ*ϵ*λt*ρ - βa*γ*(ϵ + θt - θt*μt)*τa +δ*θa*(ϵ + θt - θt*μt)*τa) + βa*δ*ϵ*(ϵ + θt)*ρ*ψt)) +(lt*αa*βa*δ*ϵ*(ϵ + θt)*ρ - αa*αt*βa*δ*ϵ*λt*ρ + αa*αt*βa*γ*ϵ*τa - αa*αt*δ*ϵ*θa*τa + αa*αt*βa*γ*θt*τa - αa*αt*δ*θa*θt*τa - αa*αt*δ*ϵ*θa*μa*τa - αa*αt*δ*θa*θt*μa*τa - αa*αt*βa*γ*θt*μt*τa + αa*αt*δ*θa*θt*μt*τa + αa*αt*δ*θa*θt*μa*μt*τa + αt*βa*δ*ϵ*θa*μa*ρ*ψa + αa*βa*δ*(ϵ**2)*ρ*ψt + αa*βa*δ*ϵ*θt*ρ*ψt + e*(la*αt*βa*(βa*γ - δ*θa)*(ϵ + θt - θt*μt) -αt*βa*(βa*γ*θt*(-1 + μt) + δ*θa*(ϵ + θt - θt*μt))*ψa -αa*(ϵ + θt - θt*μt)*(αt*(βa*δ*λa + βa*γ*ηa*τa - 2*δ*ηa*θa*τa) + βa*δ*(ϵ + θt)*ψt))**2)))/(e*αa*αt*βa*(ϵ + θt - θt*μt)*(e*βa*γ - (βa*γ + δ*θa*(-1 + μa))*ρ)),
                                                       Na = -0.5*(-(e*la*αt*(βa**2)*γ*ϵ) + e*la*αt*βa*δ*ϵ*θa - e*la*αt*(βa**2)*γ*θt + e*la*αt*βa*δ*θa*θt + e*αa*αt*βa*δ*ϵ*λa + e*αa*αt*βa*δ*θt*λa + e*la*αt*(βa**2)*γ*θt*μt - e*la*αt*βa*δ*θa*θt*μt - e*αa*αt*βa*δ*θt*λa*μt - lt*αa*βa*δ*(ϵ**2)*ρ - lt*αa*βa*δ*ϵ*θt*ρ + αa*αt*βa*δ*ϵ*λt*ρ - αa*αt*βa*γ*ϵ*τa + e*αa*αt*βa*γ*ϵ*ηa*τa + αa*αt*δ*ϵ*θa*τa - 2*e*αa*αt*δ*ϵ*ηa*θa*τa - αa*αt*βa*γ*θt*τa + e*αa*αt*βa*γ*ηa*θt*τa + αa*αt*δ*θa*θt*τa - 2*e*αa*αt*δ*ηa*θa*θt*τa + αa*αt*δ*ϵ*θa*μa*τa + αa*αt*δ*θa*θt*μa*τa + αa*αt*βa*γ*θt*μt*τa - e*αa*αt*βa*γ*ηa*θt*μt*τa - αa*αt*δ*θa*θt*μt*τa + 2*e*αa*αt*δ*ηa*θa*θt*μt*τa - αa*αt*δ*θa*θt*μa*μt*τa + e*αt*βa*δ*ϵ*θa*ψa - e*αt*(βa**2)*γ*θt*ψa + e*αt*βa*δ*θa*θt*ψa + e*αt*(βa**2)*γ*θt*μt*ψa - e*αt*βa*δ*θa*θt*μt*ψa - αt*βa*δ*ϵ*θa*μa*ρ*ψa + e*αa*βa*δ*(ϵ**2)*ψt + 2*e*αa*βa*δ*ϵ*θt*ψt + e*αa*βa*δ*(θt**2)*ψt - e*αa*βa*δ*ϵ*θt*μt*ψt - e*αa*βa*δ*(θt**2)*μt*ψt - αa*βa*δ*(ϵ**2)*ρ*ψt - αa*βa*δ*ϵ*θt*ρ*ψt + sqrt(4*αa*αt*δ*(αa*(ϵ + θt - θt*μt)*τa - βa*ϵ*ρ*ψa + e*(ϵ + θt - θt*μt)*(la*βa - αa*ηa*τa + βa*ψa))*(-(e*(lt*(βa**2)*γ*ϵ*(ϵ + θt) - αt*((βa**2)*γ*(ϵ*(λa + λt) - θt*λa*(-1 + μt)) +δ*ηa*(θa**2)*(ϵ + θt - θt*μt)*τa -βa*θa*(ϵ + θt - θt*μt)*(δ*λa + γ*ηa*τa)) + βa*(ϵ + θt)*(βa*γ*θt*(-1 + μt) + δ*θa*(ϵ + θt - θt*μt))*ψt)) + θa*μa*(lt*βa*δ*ϵ*(ϵ + θt)*ρ -αt*(βa*δ*ϵ*λt*ρ - βa*γ*(ϵ + θt - θt*μt)*τa +δ*θa*(ϵ + θt - θt*μt)*τa) + βa*δ*ϵ*(ϵ + θt)*ρ*ψt)) +(lt*αa*βa*δ*ϵ*(ϵ + θt)*ρ - αa*αt*βa*δ*ϵ*λt*ρ + αa*αt*βa*γ*ϵ*τa -αa*αt*δ*ϵ*θa*τa + αa*αt*βa*γ*θt*τa - αa*αt*δ*θa*θt*τa -αa*αt*δ*ϵ*θa*μa*τa - αa*αt*δ*θa*θt*μa*τa - αa*αt*βa*γ*θt*μt*τa +αa*αt*δ*θa*θt*μt*τa + αa*αt*δ*θa*θt*μa*μt*τa + αt*βa*δ*ϵ*θa*μa*ρ*ψa +αa*βa*δ*(ϵ**2)*ρ*ψt + αa*βa*δ*ϵ*θt*ρ*ψt +e*(la*αt*βa*(βa*γ - δ*θa)*(ϵ + θt - θt*μt) -αt*βa*(βa*γ*θt*(-1 + μt) + δ*θa*(ϵ + θt - θt*μt))*ψa -αa*(ϵ + θt - θt*μt)*(αt*(βa*δ*λa + βa*γ*ηa*τa - 2*δ*ηa*θa*τa) + βa*δ*(ϵ + θt)*ψt))**2)))/(αa*αt*δ*(-(αa*(ϵ + θt - θt*μt)*τa) + βa*ϵ*ρ*ψa -e*(ϵ + θt - θt*μt)*(la*βa - αa*ηa*τa + βa*ψa))),
                                                       ha_prod = βa*e*Ha*Pa+δ*ρ*Lt*Ha,
                                                       pt_prod = αt*Pt*Nt,
                                                       pa_prod = αa*Pa*Na,
                                                       pt_flux = ϵ*Pt,
                                                       nt_flux = ψt*Nt,
                                                       na_flux = ψa*Na)
  return(sensitivity_params)
  
}

# function to select only feasible and stable parameter combinations
feasibility_check <- function(parameter_df){
  
  # remove unfeasible equilibrium values
  feasible_params <- parameter_df %>% 
    filter(Pt > 0 & 
             Nt > 0 & 
             Lt > 0 & 
             Ha > 0 & 
             Pa > 0 & 
             Na > 0)
  
  return(feasible_params)
}

stability_check <- function(feasible_params){
  # remove parameter combinations that result in unstable equilibria
  stability_df <- NULL
  
  for (i in 1:nrow(feasible_params)) {
    # fetch parameter values from the equilibrium simulations df
    αt = feasible_params$αt[i]
    αa = feasible_params$αa[i]
    βa = feasible_params$βa[i]
    e = feasible_params$e[i]
    δ = feasible_params$δ[i]
    γ = feasible_params$γ[i]
    ϵ = feasible_params$ϵ[i]
    ηa = feasible_params$ηa[i]
    θt = feasible_params$θt[i]
    θa = feasible_params$θa[i]
    λt = feasible_params$λt[i]
    λa = feasible_params$λa[i]
    μt = feasible_params$μt[i]
    μa = feasible_params$μa[i]
    ρ = feasible_params$ρ[i]
    τa = feasible_params$τa[i]
    ψt = feasible_params$ψt[i]
    ψa = feasible_params$ψa[i]
    lt = feasible_params$lt[i]
    la = feasible_params$la[i]
    
    # Put in the solutions to the equilibrium values here:
    dPt = feasible_params$Pt[i]
    dNt = feasible_params$Nt[i]
    dLt = feasible_params$Lt[i]
    dHa = feasible_params$Ha[i]
    dPa = feasible_params$Pa[i]
    dNa = feasible_params$Na[i]
    
    # create the Jacobian from Mathematica
    jacob <- rbind(c(dNt*αt-ϵ-θt, dPt*αt, 0, 0, 0, 0), 
                   c(-(dNt*αt)+θt*μt, -lt-dPt*αt-ψt, 0, 0, 0, ψa), 
                   c(ϵ, 0, -γ-dHa*δ, -(dLt*δ), 0, 0), 
                   c(0, 0, dHa*δ*ρ, e*dPa*βa+dLt*δ*ρ-τa, e*dHa*βa, 0), 
                   c(0, 0, 0, -(dPa*βa), dNa*αa-dHa*βa-θa, dPa*αa), 
                   c(0, ψt, γ, ηa*τa, -(dNa*αa)+θa*μa, -la-dPa*αa-ψa))
    
    
    # Check to see if the real part of the leading eigenvalue is less than 0 or not - if it is than the system is stable.
    if (max(Re(base::eigen(jacob)$values)) < 0 ){ 
      stability <- "stable"
    } else {
      stability <- "unstable"
    }
    
    stability_df <- rbind(stability_df, data.frame(αt,
                                                    αa, 
                                                    βa, 
                                                    e,
                                                    δ, 
                                                    γ,
                                                    ϵ,
                                                    ηa, 
                                                    θt, 
                                                    θa,
                                                    λt, 
                                                    λa, 
                                                    μt,
                                                    μa, 
                                                    ρ,
                                                    τa, 
                                                    ψt,
                                                    ψa,
                                                    lt,
                                                    la,
                                                    dPt, 
                                                    dNt, 
                                                    dLt,
                                                    dHa, 
                                                    dPa, 
                                                    dNa, 
                                                    EiV1 = Re(eigen(jacob)$values[1]),
                                                    EiV2 = Re(eigen(jacob)$values[2]),
                                                    EiV3 = Re(eigen(jacob)$values[3]),
                                                    EiV4 = Re(eigen(jacob)$values[4]),
                                                    EiV5 = Re(eigen(jacob)$values[5]),
                                                    EiV6 = Re(eigen(jacob)$values[6]),
                                                    maxEv = max(Re(base::eigen(jacob)$values)),
                                                    stability = stability))
  }
  
  # join stability and feasible_params together and filter out rows that produce unstable or negative equilibrium values
  stable_params <- feasible_params %>% 
    left_join(stability_df %>% select(-dPt,
                                      -dNt,
                                      -dLt,
                                      -dHa,
                                      -dPa,
                                      -dNa,
                                      -EiV1,
                                      -EiV2,
                                      -EiV3,
                                      -EiV4,
                                      -EiV5,
                                      -EiV6,
                                      -maxEv), by = c("αt",
                                                      "αa",
                                                      "βa",
                                                      "e",
                                                      "δ",
                                                      "γ",
                                                      "ϵ",
                                                      "ηa",
                                                      "θt",
                                                      "θa",
                                                      "λt",
                                                      "λa",
                                                      "μt",
                                                      "μa",
                                                      "ρ",
                                                      "τa",
                                                      "ψt",
                                                      "ψa",
                                                      "lt",
                                                      "la")) %>% 
    filter(stability == "stable") %>%
    select(-stability)
  
  return(stable_params)
}

tic("forestry disturbance simulation")
# 1) Forestry/insect outbreak + moose
sensitivity_params_θt_μt <- simulate_disturbance(master_df=(param_simulations_stable %>% head(10000)),
                                             param_number=9,  # θt 
                                             min_param_val=0,
                                             max_param_val=10,
                                             increase_val=0.5,
                                             expand_length=21,
                                             param_number2=13,   # μt
                                             min_param_val2=0,
                                             max_param_val2=1, 
                                             increase_val2=0.05,
                                             expand_length2=21)

feasible_data_θt_μt <- feasibility_check(sensitivity_params_θt_μt)
stable_data_θt_μt <- stability_check(feasible_data_θt_μt)
write.csv(stable_data_θt_μt, "output/stable_data_θt_μt_10000.csv", row.names=FALSE)
toc()

tic("generate data for forest simulation plot")
# transform data to matrix form
surface_data_forestry <- stable_data_θt_μt %>% 
  select(θt, μt, Ha, Pa, ha_prod, pa_prod, Na, Nt) %>% 
  group_by(θt, μt) %>% 
  summarise_at(vars(Ha, Pa, ha_prod, pa_prod, Na, Nt),
               funs(mean(., na.rm=TRUE))) %>% 
  ungroup()

write.csv(surface_data_forestry, "output/surface_data_forestry_10000.csv", row.names=FALSE)

p1 <- ggplot(surface_data_forestry, aes(θt, μt, z=Ha)) +
  geom_contour_filled() +
  theme_classic() + 
   guides(fill=guide_legend(title="Ha density")) + 
  xlab("") + 
  ylab("proportion recycled (μt)")

p2 <- ggplot(surface_data_forestry, aes(θt, μt, z=ha_prod)) +
  geom_contour_filled() +
  theme_classic() + 
  guides(fill=guide_legend(title="Ha productivity")) + 
  xlab("") + 
  ylab("")

p3 <- ggplot(surface_data_forestry, aes(θt, μt, z=Pa)) +
  geom_contour_filled() +
  theme_classic() + 
  guides(fill=guide_legend(title="Pa density")) + 
  xlab("death rate of Pt (θt)") + 
  ylab("proportion recycled (μt)")

p4 <- ggplot(surface_data_forestry, aes(θt, μt, z=pa_prod)) +
  geom_contour_filled() +
  theme_classic() + 
  guides(fill=guide_legend(title="Pa productivity")) + 
  xlab("death rate of Pt (θt)") + 
  ylab("")

forestry_surface_plots <- grid.arrange(p1, p2, p3, p4, nrow=2)
toc()
ggsave(file = "output/surafce_plots_forestry_10000.svg", plot=forestry_surface_plots, width=12, height=10)

p5 <- ggplot(surface_data_forestry, aes(θt, μt, z=Na)) +
  geom_contour_filled() +
  theme_classic() + 
  guides(fill=guide_legend(title="Na density")) + 
  xlab("death rate of Pt (θt)") + 
  ylab("proportion recycled (1-μt)") +
  scale_y_reverse()

p6 <- ggplot(surface_data_forestry, aes(θt, μt, z=Nt)) +
  geom_contour_filled() +
  theme_classic() + 
  guides(fill=guide_legend(title="Nt density")) + 
  xlab("death rate of Pt (θt)") + 
  ylab("") +
  scale_y_reverse()
forestry_nutrient_surface_plots <- grid.arrange(p5, p6, nrow=1)
ggsave(file = "output/surafce_plots_nutrients_forestry_10000.svg", plot=forestry_nutrient_surface_plots, width=12, height=6)

tic("atv trail simulation")
# 2) ATV trails
sensitivity_params_βa_αa <- simulate_disturbance(master_df=(param_simulations_stable %>% head(10000)),
                                             param_number=2,  # αa
                                             min_param_val=0,
                                             max_param_val=10,
                                             increase_val=0.5,
                                             expand_length=21,
                                             param_number2=3,   # βa
                                             min_param_val2=0,
                                             max_param_val2=10, 
                                             increase_val2=0.5,
                                             expand_length2=21)

feasible_data_βa_αa <- feasibility_check(sensitivity_params_βa_αa)
stable_data_βa_αa <- stability_check(feasible_data_βa_αa)
write.csv(stable_data_βa_αa, "output/stable_data_βa_αa_10000.csv", row.names=FALSE)
toc()

tic("plot atv simulation")
# transform data to matrix form
surface_data_atv <- stable_data_βa_αa %>% 
  select(βa, αa, Ha, Pa, ha_prod, pa_prod, Na, Nt) %>% 
  group_by(βa, αa) %>% 
  summarise_at(vars(Ha, Pa, ha_prod, pa_prod, Na, Nt),
               funs(mean(., na.rm=TRUE))) %>% 
  ungroup()

write.csv(surface_data_atv, "output/surface_data_atv_10000.csv", row.names=FALSE)

p1 <- ggplot(surface_data_atv, aes(βa, αa, z=Ha)) +
  geom_contour_filled() +
  theme_classic() + 
  guides(fill=guide_legend(title="Ha density")) + 
  xlab("") + 
  ylab("Pa uptake (αa)") + 
  scale_x_reverse() + 
  scale_y_reverse()

p2 <- ggplot(surface_data_atv, aes(βa, αa, z=ha_prod)) +
  geom_contour_filled() +
  theme_classic() + 
  guides(fill=guide_legend(title="Ha productivity")) + 
  xlab("") + 
  ylab("") + 
  scale_x_reverse() + 
  scale_y_reverse()

p3 <- ggplot(surface_data_atv, aes(βa, αa, z=Pa)) +
  geom_contour_filled() +
  theme_classic() + 
  guides(fill=guide_legend(title="Pa density")) +
  xlab("Ha uptake (βa)") + 
  ylab("Pa uptake (αa)") + 
  scale_x_reverse() + 
  scale_y_reverse()

p4 <- ggplot(surface_data_atv, aes(βa, αa, z=pa_prod)) +
  geom_contour_filled() +
  theme_classic() + 
  guides(fill=guide_legend(title="Pa productivity")) + 
  xlab("Ha uptake (βa)") + 
  ylab("") + 
  scale_x_reverse() + 
  scale_y_reverse()

atv_surface_plots <- grid.arrange(p1, p2, p3, p4, nrow=2)
ggsave(file = "output/surafce_plots_atv_10000.svg", plot=atv_surface_plots, width=12, height=10)
toc() 

p5 <- ggplot(surface_data_atv, aes(βa, αa, z=Na)) +
  geom_contour_filled() +
  theme_classic() + 
  guides(fill=guide_legend(title="Na density")) + 
  xlab("Ha uptake (βa)") + 
  ylab("Pa uptake (αa)") + 
  scale_x_reverse() + 
  scale_y_reverse()

p6 <- ggplot(surface_data_atv, aes(βa, αa, z=Nt)) +
  geom_contour_filled() +
  theme_classic() + 
  guides(fill=guide_legend(title="Nt density")) +
  xlab("Decreasing Ha uptake (1-βa)") + 
  ylab("") + 
  scale_x_reverse() + 
  scale_y_reverse()

atv_nutrient_surface_plots <- grid.arrange(p5, p6, nrow=1)
ggsave(file = "output/surafce_plots_atv_nutrients_10000.svg", plot=atv_nutrient_surface_plots, width=12, height=6)

