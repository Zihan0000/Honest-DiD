rm(list = ls())
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
source(here("C:/Users/22133/Desktop/master thesis/HonestParallelTrends_RESTUD_ReplicationFiles/HonestParallelTrends_RESTUD_ReplicationFiles/Code/simulationStudy/DeltaSD_excessLength_functions.R"))
set.seed(10009)
####DGP
numPrePeriods<-5
numPostPeriods<-10
p<-numPrePeriods+numPostPeriods+1
N=200
grid  <- ffscb::make_grid(p, rangevals=c(-1,1))
eps <- matrix(rnorm(p * N,0,1), nrow = p, ncol = N)#6
phi_t <- meanf_shift((grid+1)/2, 0) * 1#2.5
lambda_i <- rnorm(N, mean=0, sd=0.5)#0.5
D_i <- sample(0:1, size=N, replace=TRUE)
one <- diag(1,nrow = p)
beta<-rep(0,p)
#beta<-c(-0.2,-2,-0.3,-3,-3,-6,-5,-11,-12,-14,-17,-19,-18,-20,-20,-18)
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


##define the M and L for SD
Mvec = as.vector(seq(from = 0, to = 5, by = 1))
l_vec = rep(1, numPostPeriods)/numPostPeriods
truebeta<-rep(0,p-1)

#run simulation
run_simulation_sd <- function(num_simulations) {
  betahat_m= rmvnorm(num_simulations, mean = truebeta, sigma = sigma)
  Result_PT_sd_Ave<-list()
  # Run the simulation for each iteration
  for (i in 1:num_simulations) {
    betahat=betahat_m[i,]
    Result_PT_sd_Ave[[i]]<-compareCI_EL_DeltaSD_sim(1,trueBeta = truebeta,
                                                    betahat=betahat,
                                                    numPrePeriods = numPrePeriods,
                                                    numPostPeriods =numPostPeriods,
                                                    sigma=sigma,
                                                    Mvec = Mvec,
                                                    l_vec=l_vec,
                                                    alpha=0.05,
                                                    rescalePermanently = T)
  }
  #Combine results from all simulations into a single data frame
  combined_data <- bind_rows(Result_PT_sd_Ave,.id="i")
  # Summarize the mean 'EL' for each combination of 'M' and 'CItype' across all simulations
  result_sd_EL <- combined_data %>%
    group_by(M, CItype) %>%
    summarize(mean_EL = mean(EL),sd=sd(EL))
  return(result_sd_EL)
}

example<-run_simulation_sd(100)
library(furrr)
##then calculate the optimal EL:
optimalEL_df<- computeOptimalEL_DeltaSD_Mvec(trueBeta = truebeta, 
                                             sigma = sigma, 
                                             Mvec = Mvec, 
                                             l_vec = l_vec,
                                             numPrePeriods = numPrePeriods, 
                                             numPostPeriods = numPostPeriods, 
                                             alpha = 0.05,
                                             rescalePermanently = T)
colnames(optimalEL_df)[2] <- "Opt_EL"


#ratio
optimalEL_df_selecte <- optimalEL_df %>%
  dplyr::select(M, Opt_EL)
merged_SD<- merge(example, optimalEL_df, by = "M", all.x = TRUE)
merged_SD <- merged_SD %>%
  mutate(efficiency_sd = Opt_EL / mean_EL)
saveRDS(merged_SD, file = "merged_SD_Ave.rds")

# Construct simulation results -------------------------------------------------------
# Conditional, FLCI, C-LF Hybrid
library(ggplot2)
source(here("C:/Users/22133/Desktop/master thesis/HonestParallelTrends_RESTUD_ReplicationFiles/HonestParallelTrends_RESTUD_ReplicationFiles/Code/UtilityFunctions/fte_theme.R"))
plot = ggplot(merged_SD %>% filter(CItype.x != "C-F Hybrid"), 
              aes(x = M, y = efficiency_sd)) +
  geom_line(aes(colour = CItype.x)) + geom_point(aes(colour = CItype.x, shape = CItype.x)) +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25), limits = c(0, 1.05)) +
  scale_colour_manual(values = c("#377EB8", "#4DAF4A", "#984EA3")) +
  fte_theme() + 
  theme(legend.title = element_blank(), legend.position = "none",
        plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle(latex2exp::TeX("$\\Delta = \\Delta^{SD}(M)$, $\\theta = \\bar{\\tau}$")) +
  labs(x = latex2exp::TeX("$M$"), y = "Excess length ratio") 

plot


ggsave(plot, 
       file = here("Figures/sd_avg_withtag.png"), 
       width = 5, height = 3.5, units = 'in')

# Clear environment
rm(list = ls())








