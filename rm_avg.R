# Description:
#   This code implements ratio of excess length(EL) for the set Delta^{RM}(Mbar).
#   The script is organized into the following parts:
#     1) Data Generating Process -- Generate data
#     2) Sigma estimation by TWFE -- use TWFE to estimate sigma and betahat
#     3) Simulation helper functions -- implement functions to conduct the excess length simulations.
#     4) Optimal length functions -- implement helper functions to compute the bound on the optimal length.
#     5) Ratio of EL -- the ratio between excess length and Optimal length.

rm(list = ls())
# Libraries
library(ROI)
library(foreach)
library(doParallel)
library(here)
library(CVXR)
library("purrr")
library(dplyr)
library("ffscb")
library(HonestDiD)
library(mvtnorm)

#source from the original code of Prof. Rambachan and Roth
source(here("C:/Users/22133/Desktop/master thesis/HonestParallelTrends_RESTUD_ReplicationFiles/HonestParallelTrends_RESTUD_ReplicationFiles/Code/simulationStudy/DeltaRM_excessLength_functions.R"))
set.seed(10080)

####DGP
numPrePeriods<-3
numPostPeriods<-5
p<-numPrePeriods+numPostPeriods+1
N=200
grid  <- ffscb::make_grid(p, rangevals=c(-1,1))
eps <- matrix(rnorm(p * N,0,1), nrow = p, ncol = N)
phi_t <- meanf_shift((grid+1)/2, 0) * 1
lambda_i <- rnorm(N, mean=0, sd=0.5)
D_i <- sample(0:1, size=N, replace=TRUE)
one <- diag(1,nrow = p)
beta<-c(rep(0,numPrePeriods-1),0.5,rep(0,numPostPeriods+1))
Y <-  t(lambda_i + t(phi_t + beta * (outer(diag(one), D_i)) + eps)) 
data <- data.frame(
  individual_id = rep(1:N, each = p),
  year = rep(0:(p-1), times = N),
  treatment = rep(D_i, each = p),
  Y_it = as.vector(Y)
)

###TWFE
twfe_results <- fixest::feols(Y_it ~ i(year, treatment, ref = 3)|individual_id + year, data=data)
betahat <- summary(twfe_results)$coefficients
sigma <- summary(twfe_results)$cov.scaled

#debug
dataex<-surveyOfPapers_resultsObjectList[[5]]
sigma<-dataex$sigma
Mbarvec<-c(1)
l_vec<-c(1,rep(0,3))

##define the Mbar and L for RM
Mbarvec = c(1)
l_vec = rep(1, numPostPeriods)/numPostPeriods

#simulation part for the excess lenghth of RM

run_simulation_rm<- function(num_simulations) {
  EL_df_results = tibble()
  # Run the simulation for each iteration
  for (multiplier in c(0, 1, 2, 3)) {
    # Create true beta and extract simulated data associated with true beta.
    trueBeta = multiplier * basisVector(index = numPrePeriods, numPrePeriods + numPostPeriods)
    betahat_m= rmvnorm(num_simulations, mean = trueBeta, sigma = sigma)
    
    for (i in 1:num_simulations) {
      betahat=betahat_m[i,]
      temp<-compareCI_EL_DeltaRM_sim(1,trueBeta = trueBeta,
                                                     betahat=betahat,
                                                     numPrePeriods = numPrePeriods,
                                                     numPostPeriods =numPostPeriods,
                                                     sigma=sigma,
                                                     Mbarvec = Mbarvec,
                                                     l_vec=l_vec,
                                                     alpha=0.05,
                                                     rescalePermanently = T,
                                                     gridPoints = 1000, grid.ub = NA, grid.lb = NA)
      temp$Mbarvec = Mbarvec
      temp$multiplier = multiplier
      EL_df_results <- bind_rows(EL_df_results, temp)
    }
  }
  result_rm_EL <- EL_df_results %>%
    group_by(CItype,multiplier) %>%
    summarize(mean_EL = mean(EL), sd = sd(EL))
  return(result_rm_EL)
  
}
example_rm<-run_simulation_rm(100)
colnames(example_rm)[2] <- "pulse_multiplier"
saveRDS(example_rm, file = "example_rm_avg.rds")
#calculate the optimal excess length for RM
optEL_df_results = tibble()
for (multiplier in c(0, 1, 2, 3)) {

  trueBeta = multiplier * basisVector(index = numPrePeriods, numPrePeriods + numPostPeriods)
  
  # Compute excess length for optimal CI
  optimalEL_df = computeOptimalEL_DeltaRM_Mbarvec(trueBeta = trueBeta, sigma = sigma, 
                                                  Mbarvec = Mbarvec, l_vec = l_vec,
                                                  numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods, 
                                                  alpha = 0.05, rescalePermanently = T)
  optimalEL_df = optimalEL_df %>% 
    group_by(Mbar) %>%
    mutate(
      pulse_multiplier = rep(multiplier, 1)
    )
  
  # Add results to summary.
  optEL_df_results <- bind_rows(
    optEL_df_results, optimalEL_df
  )
}
colnames(optEL_df_results)[2] <- "Opt_EL"
optEL_df_results
saveRDS(optEL_df_results, file = "optEL_df_results_avg.rds")
#ratio
merged_RM<- merge(example_rm, optEL_df_results, by = "pulse_multiplier", all.x = TRUE)
merged_RM <- merged_RM %>%
  mutate(efficiency_rm = Opt_EL / mean_EL) 
saveRDS(merged_RM, file = "merged_RM_avg.rds")
## Construct simulation results -------------------------------------------------------
# Conditional, C-LF Hybrid
# Construct Simulation Results for Mbar = 1
library(ggplot2)
source(here("C:/Users/22133/Desktop/master thesis/HonestParallelTrends_RESTUD_ReplicationFiles/HonestParallelTrends_RESTUD_ReplicationFiles/Code/UtilityFunctions/fte_theme.R"))
plot = ggplot(merged_RM %>% filter(Mbar == 1), 
              aes(x = pulse_multiplier, y = efficiency_rm, color = CItype.x)) +
  geom_line() + geom_point(aes(shape = CItype.x)) +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25), limits = c(0, 1.05)) +
  fte_theme() + scale_colour_manual(values = c("#377EB8", "#4DAF4A")) +
  scale_shape_manual(values=c(17, 15)) + 
  theme(legend.title = element_blank(), 
        legend.position = c(0.8, 0.2),
        legend.text = element_text(size = 7),
        legend.background = element_rect(colour = "transparent", fill = alpha("red", 0)),
        plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle(latex2exp::TeX("$\\Delta^{RM}(\\bar{M})$, $\\theta = \\bar{\\tau}$")) +
  labs(x = latex2exp::TeX("$\\delta_{-1}$"), y = "Excess length ratio", 
       colour = "", shape = "", linetype = "")
plot
ggsave(plot, 
       file = here("Figures/Figure_rm_avg.png"), 
       width = 5, height = 3.5, units = 'in')
