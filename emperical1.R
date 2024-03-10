library(did)
library(here)
library(dplyr)
library(did)
library(haven)
library(ggplot2)
library(fixest)
library(HonestDiD)
data(mpdta)#the data is from callway's paper.
head(mpdta)
data <- subset(mpdta, first.treat == 2004 | first.treat == 0)
pre_period<-1
post_period=3

twfe_results <- fixest::feols(lemp ~ i(year, treat, ref = 2004) |countyreal+year, 
                              cluster = "countyreal",
                              data = data)
betahat <- summary(twfe_results)$coefficients #save the coefficients
sigma <- summary(twfe_results)$cov.scaled #save the covariance matrix
fixest::iplot(twfe_results,pt.col= "#C12928",col = "#27c2b2",
              group.par = list(lwd = 1, line = 3, tcl = 0.75),
              pt.lwd= 2,
              ci.lwd=2,
              zero=F,
              ref.line = F,
              ci.width=0.15,
              main = "Effect on teen employment(log)",
              lwd = 1.5)

#RM results:
delta_rm_results <- 
  HonestDiD::createSensitivityResults_relativeMagnitudes(
    betahat = betahat, #coefficients
    sigma = sigma, #covariance matrix
    numPrePeriods = 1, #num. of pre-treatment coefs
    numPostPeriods = 3, #num. of post-treatment coefs
    Mbarvec = seq(0,2,by=0.5)#values of Mbar
  )
#ols
originalResults <- HonestDiD::constructOriginalCS(betahat = betahat,
                                                  sigma = sigma,
                                                  numPrePeriods = 1,
                                                  numPostPeriods = 3)
#compare ols and RM results
HonestDiD::createSensitivityPlot_relativeMagnitudes(delta_rm_results, originalResults)
#SD
delta_sd_results <- 
  HonestDiD::createSensitivityResults(betahat = betahat,
                                      sigma = sigma,
                                      numPrePeriods = 1,
                                      numPostPeriods = 3,
                                      Mvec = seq(from = 0, to = 0.05, by =0.01))
#compare
createSensitivityPlot(delta_sd_results, originalResults)
#talk about the situation 
#####for lvec=average
delta_rm_results_avg <- 
  HonestDiD::createSensitivityResults_relativeMagnitudes(betahat = betahat,
                                                         sigma = sigma,
                                                         numPrePeriods = 1,
                                                         numPostPeriods = 3, Mbarvec = seq(0,2,by=0.5),
                                                         l_vec = c(1/3,1/3,1/3))

originalResults_avg <- HonestDiD::constructOriginalCS(betahat = betahat,
                                                      sigma = sigma,
                                                      numPrePeriods = 1,
                                                      numPostPeriods = 3,
                                                      l_vec = c(1/3,1/3,1/3))
HonestDiD::createSensitivityPlot_relativeMagnitudes(delta_rm_results_avg, originalResults_avg)
##for SD
delta_sd_results_avg <- 
  HonestDiD::createSensitivityResults(betahat = betahat,
                                      sigma = sigma,
                                      numPrePeriods = 1,
                                      numPostPeriods = 3,
                                      Mvec = seq(from = 0, to = 0.05, by =0.01),
                                      l_vec = c(1/3,1/3,1/3))
#compare
createSensitivityPlot(delta_sd_results_avg, originalResults_avg)
