# Description:
#   This code implements ratio of excess length(EL) for the set Delta^{SDPB}(M).
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
source(("C:/Users/22133/Desktop/master thesis/HonestParallelTrends_RESTUD_ReplicationFiles/HonestParallelTrends_RESTUD_ReplicationFiles/Code/simulationStudy/DeltaSDPB_excessLength_functions.R"))
set.seed(1010)
####DGP
numPrePeriods<-5
numPostPeriods<-10
p<-numPrePeriods+numPostPeriods+1
N=200
grid  <- ffscb::make_grid(p, rangevals=c(-1,1))
eps <- matrix(rnorm(p * N,0,1), nrow = p, ncol = N)
phi_t <- meanf_shift((grid+1)/2, 0)*1
lambda_i <- rnorm(N, mean=0, sd=0.5)
D_i <- sample(0:1, size=N, replace=TRUE)
one <- diag(1,nrow = p)
#beta<-c(-0.03,-0.01,-0.02,-0.01,-0.05,0,0.02,0.025,0.036,0.01,0.03,0.06,0.07,0.1,0.2,0.23)
beta<-c(rep(0,p))
#beta<-c(0.2,2,0.3,3,3,6,5,11,12,14,17,19,18,20,20,18)
Y <-  t(lambda_i + t(phi_t + beta * (outer(diag(one), D_i)) + eps)) 
data <- data.frame(
  individual_id = rep(1:N, each = p),
  year = rep(0:(p-1), times = N),
  treatment = rep(D_i, each = p),
  Y_it = as.vector(Y)
)

###TWFE
twfe_results <- fixest::feols(Y_it ~ i(year, treatment, ref = 5)|individual_id + year, data=data)
betahat <- summary(twfe_results)$coefficients
sigma <- summary(twfe_results)$cov.scaled


##define the Mbar and L for SDPB
Mvec = as.vector(seq(from = 0, to = 5, by = 1))
l_vec=c(1,rep(0,9))
l_vec = rep(1, numPostPeriods)/numPostPeriods
truebeta<-rep(0,p-1)
#simulation part for the excess lenghth of SDPB
run_simulation_SDPB<- function(num_simulations) {
  betahat_m= rmvnorm(num_simulations, mean = truebeta, sigma = sigma)
  Result_PT_sdpb_Ave<-list()
  # Run the simulation for each iteration
  for (i in 1:num_simulations) {
    betahat=betahat_m[i,]
    Result_PT_sdpb_Ave[[i]]<-compareCI_EL_DeltaSDPB_sim(1,trueBeta = truebeta,
                                                    betahat=betahat,
                                                    numPrePeriods = numPrePeriods,
                                                    numPostPeriods =numPostPeriods,
                                                    sigma=sigma,
                                                    Mvec = Mvec,
                                                    l_vec=l_vec,
                                                    alpha=0.05,
                                                    rescalePermanently = T
                                                    )
  }
  #Combine results from all simulations into a single data frame
  combined_data <- bind_rows(Result_PT_sdpb_Ave,.id="i")
  #Summarize the mean 'EL' for each combination of 'M' and 'CItype' across all simulations
  result_sdpb_EL <- combined_data %>%
    group_by(M, CItype) %>%
    summarize(mean_EL = mean(EL),sd=sd(EL))
  return(result_sdpb_EL)
}
example_sdpb<-run_simulation_SDPB(100)
#calculate the optimal excess length for SDPB
library(furrr)
optimal_EL_sdpb<-computeOptimalEL_DeltaSDPB_Mvec(trueBeta=truebeta, sigma=sigma, Mvec=Mvec, l_vec=l_vec,
                                                numPrePeriods=numPrePeriods, numPostPeriods=numPostPeriods, 
                                                alpha=0.05, rescalePermanently=T)

colnames(optimal_EL_sdpb)[2] <- "Opt_EL"
#saveRDS(optimal_EL_sdpb, file = "optimalEL_df_SD_Ave.rds")
loaded_data <- readRDS("merged_SDPB_1.rds.rds")


#ratio
#optimalEL_sdpb_selecte <- optimal_EL_sdpb %>%
  #dplyr::select(M, Opt_EL)
merged_SDPB<- merge(example_sdpb, optimal_EL_sdpb, by = "M", all.x = TRUE)
merged_SDPB <- merged_SDPB %>%
  mutate(efficiency_sdpb = Opt_EL / mean_EL)
saveRDS(merged_SDPB, file = "merged_SDPB_1.rds")
readRDS("merged_SDPB_1.rds")
# Construct simulation results -------------------------------------------------------
merged_SDPB <- merged_SDPB_1 %>%
  mutate(CItype_renamed = case_when(
    CItype.x == "C-F w/ pos. bias" ~ "C-F Hybrid",
    CItype.x == "C-LF w/ pos. bias" ~ "C-LF Hybrid",
    CItype.x == "Conditional w/ pos. bias" ~ "Conditional",
    CItype.x == "FLCI" ~ "FLCI"
  ))
library(ggplot2)
source(here("Code/UtilityFunctions/fte_theme.R"))
# Conditional, FLCI, C-LF Hybrid
plot = ggplot(merged_SDPB %>% filter(CItype_renamed != "C-F Hybrid"), 
              aes(x = M, y = efficiency_sdpb)) +
  geom_line(aes(colour = CItype_renamed)) + geom_point(aes(colour = CItype_renamed, shape = CItype_renamed)) +
  scale_y_continuous(breaks = seq(from = 0, to = 1.6, by = 0.25), limits = c(0, 1.1)) +
  scale_colour_manual(values = c("#377EB8", "#4DAF4A", "#984EA3")) +
  fte_theme() + 
  theme(legend.title = element_blank(), 
        legend.position = c(0.9, 0.65),
        legend.text = element_text(size = 6),
        legend.background = element_rect(colour = "transparent", fill = alpha("red", 0)),
        plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle(latex2exp::TeX("$\\Delta = \\Delta^{SDPB}(M)$, $\\theta = {\\tau}_1$")) +
  labs(x = latex2exp::TeX("$M$"), y = " Excess length ratio")

plot
ggsave(plot, 
       file = here("Figures/sdpb_1.png"), 
       width = 5, height = 3.5, units = 'in')





