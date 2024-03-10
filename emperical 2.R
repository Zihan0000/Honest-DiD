library(did)
library(here)
library(dplyr)
library(did)
library(haven)
library(ggplot2)
library(fixest)
library(HonestDiD)
#read the list
list<-Object_list[[3]]
source(here("C:/Users/22133/Desktop/master thesis/HonestParallelTrends_RESTUD_ReplicationFiles/HonestParallelTrends_RESTUD_ReplicationFiles/Code/UtilityFunctions/applications-helper-functions.R"))
source(here("C:/Users/22133/Desktop/master thesis/HonestParallelTrends_RESTUD_ReplicationFiles/HonestParallelTrends_RESTUD_ReplicationFiles/Code/UtilityFunctions/fte_theme.R"))

eventStudyPlot <- createEventStudyPlot(list) + 
  xlab("Year") + 
  ggtitle("Per capita trans. from gov.(log)")
eventStudyPlot
ggsave(eventStudyPlot, file = here("Figures/emp2-coef.png"),
       width = 6 , height = 4)
#parameters:
betahat<-list$beta
sigma<-list$sigma
numPrePeriods<-length(list$prePeriodIndices)
numPostPeriods<-length(list$postPeriodIndices)
#RM results:
delta_rm_results <- 
  HonestDiD::createSensitivityResults_relativeMagnitudes(
    betahat = betahat, #coefficients
    sigma = sigma, #covariance matrix
    numPrePeriods = numPrePeriods, #num. of pre-treatment coefs
    numPostPeriods = numPostPeriods, #num. of post-treatment coefs
    Mbarvec = seq(0,2,by=0.5)#values of Mbar
  )
saveRDS(delta_rm_results,file = here("Output/emp2-RM.rds"))
#OLS
originalResults <- HonestDiD::constructOriginalCS(betahat = betahat,
                                                  sigma = sigma,
                                                  numPrePeriods = numPrePeriods,
                                                  numPostPeriods = numPostPeriods)
saveRDS(originalResults,file = here("Output/emp2-OLS.rds"))
#compare ols and RM results
HonestDiD::createSensitivityPlot_relativeMagnitudes(delta_rm_results, originalResults)

#SD
delta_sd_results <- 
  HonestDiD::createSensitivityResults(betahat = betahat,
                                      sigma = sigma,
                                      numPrePeriods = numPrePeriods,
                                      numPostPeriods = numPostPeriods,
                                      Mvec = seq(from = 0, to = 0.05, by =0.01))
saveRDS(delta_sd_results,file = here("Output/emp2-SD.rds"))
#compare
createSensitivityPlot(delta_sd_results, originalResults)
#talk about the situation 
#####for lvec=average
delta_rm_results_avg <- 
  HonestDiD::createSensitivityResults_relativeMagnitudes(betahat = betahat,
                                                         sigma = sigma,
                                                         numPrePeriods = numPrePeriods,
                                                         numPostPeriods = numPostPeriods, 
                                                         Mbarvec = seq(0,2,by=0.5),
                                                         l_vec = rep(1,numPostPeriods)/numPostPeriods,
                                                         grid.ub = 0.5,
                                                         grid.lb =-0.5 
                                                         ) 

saveRDS(delta_rm_results_avg,file = here("Output/emp2-RM-avg.rds"))
originalResults_avg <- HonestDiD::constructOriginalCS(betahat = betahat,
                                                      sigma = sigma,
                                                      numPrePeriods = numPrePeriods,
                                                      numPostPeriods = numPostPeriods,
                                                      l_vec = rep(1,numPostPeriods)/numPostPeriods
                                                      )
saveRDS(originalResults_avg,file = here("Output/emp2-OLS-avg.rds"))
HonestDiD::createSensitivityPlot_relativeMagnitudes(delta_rm_results_avg, originalResults_avg)
##for SD
delta_sd_results_avg <- 
  HonestDiD::createSensitivityResults(betahat = betahat,
                                      sigma = sigma,
                                      numPrePeriods = numPrePeriods,
                                      numPostPeriods = numPostPeriods,
                                      Mvec = seq(from = 0, to = 0.05, by =0.01),
                                      l_vec = rep(1,numPostPeriods)/numPostPeriods
                                      )
saveRDS(delta_sd_results_avg,file = here("Output/emp2-SD-avg.rds"))
#compare
createSensitivityPlot(delta_sd_results_avg, originalResults_avg)
#create a table to describe the absolute value

