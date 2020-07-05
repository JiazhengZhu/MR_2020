library(plotly)
library(plyr)
library(MendelianRandomization)
library(progress)

simulate = function(nsamples, alpma_max, beta, xerror, measerror, yerror){
  g = matrix(rbinom(30*nsamples,2,0.2), ncol = 30, nrow = nsamples)
  beta1 = beta[1]
  beta2 = beta[2]
  beta3 = beta[3]
  
  alpha1 = runif(30,min=0,max=alpma_max[1])
  alpha2 = runif(30,min=0,max=alpma_max[2])
  alpha3 = runif(30,min=0,max=alpma_max[3])
  #generate exposure levels x1, x2, x3 and outcome y
  x1 = rep(0,nsamples)
  x2 = rep(0,nsamples)
  x3 = rep(0,nsamples)
  y = rep(0,nsamples)
  for (i in 1:nsamples) {
    x1[i] = g[i,]%*%alpha1 + rnorm(1)*xerror[1]
    x2[i] = g[i,]%*%alpha2 + rnorm(1)*xerror[2]
    x3[i] = g[i,]%*%alpha3 + rnorm(1)*xerror[3]
    y[i] = beta1*x1[i] + beta2*x2[i] + beta3*x3[i] + rnorm(1)*yerror
  }
  
  #calculate g-outcome associations and exposure-outcome associations
  by = matrix(0,nrow = 30,ncol = 1)
  byse = matrix(0,nrow = 30,ncol = 1)
  for (i in 1:30){
    fit = lm(y~g[, i]-1)
    by[i] = coefficients(fit)
    byse[i] = coef(summary(fit))[, "Std. Error"]
  }
  
  #two-sample so regenerate exposure levels x1, x2, x3 and outcome y
  x1 = rep(0,nsamples)
  x2 = rep(0,nsamples)
  x3 = rep(0,nsamples)
  y = rep(0,nsamples)
  for (i in 1:nsamples) {
    x1[i] = g[i,]%*%alpha1 + rnorm(1)*xerror[1]
    x2[i] = g[i,]%*%alpha2 + rnorm(1)*xerror[2]
    x3[i] = g[i,]%*%alpha3 + rnorm(1)*xerror[3]
    y[i] = beta1*x1[i] + beta2*x2[i] + beta3*x3[i] + rnorm(1)*yerror
  }
  x = cbind(x1 + rnorm(1)*measerror[1],x2+ rnorm(1)*measerror[2],x3+ rnorm(1)*measerror[3])
  bx = matrix(0,nrow = 30,ncol = 3)
  bxse = matrix(0,nrow = 30,ncol = 3)
  for (i in 1:3){
    for (j in 1:30){
      fit = lm(x[,i]~g[,j]-1)
      bx[j,i] = coefficients(fit)
      bxse[j,i] = coef(summary(fit))[, "Std. Error"]
    }
  }

  #plug in the model
  mrmv = mr_mvinput(bx=bx,
                    bxse=bxse,
                    by=as.double(by),
                    byse=as.double(byse))
  
  return(list(mr_mvegger(mrmv), mr_mvivw(mrmv)))
}

egg_CI_lower = matrix(0, nrow = trials, ncol = 3)
egg_CI_upper = matrix(0, nrow = trials, ncol = 3)
egg_est = matrix(0, nrow = trials, ncol = 3)
ivw_CI_lower = matrix(0, nrow = trials, ncol = 3)
ivw_CI_upper = matrix(0, nrow = trials, ncol = 3)
ivw_est = matrix(0, nrow = trials, ncol = 3)
results = data.frame()
#true values
beta = c(0.2,-0.6,0.0)
alpha_list = list(c(1,1,1),
                  c(0.8,0.8,0.8),
                  c(0.4,0.4,0.4),
                  c(0.2,0.2,0.2))
xerror = c(0.2,0.2,0.2)
measerror_list = list(c(0,0,0),
                      c(0.1,0,0),
                      c(0.2,0,0),
                      c(0.8,0,0),
                      c(0,0.1,0),
                      c(0,0.2,0),
                      c(0,0.8,0),
                      c(0,0,0.1),
                      c(0,0,0.2),
                      c(0,0,0.8),
                      c(0.1,0.1,0.1),
                      c(0.2,0.2,0.2),
                      c(0.8,0.8,0.8))
yerror = 0.1

trials = 2000
nsamples = 1000

pb = progress_bar$new(format = " [:bar] :percent eta: :eta",total = trials*length(measerror_list)*length(alpha_list))
for (alphamax in alpha_list){
  for (measerror in measerror_list){
    egg_coverage_count = ivw_coverage_count = rep(trials,3)
    egg_pwr_count = ivw_pwr_count = rep(trials,3)
    
    for (i in 1:trials){
      pb$tick()
      m = simulate(nsamples, alphamax, beta, xerror, measerror, yerror)
      egg_CI_lower[i,] = m[[1]]@CILower.Est
      egg_CI_upper[i,] = m[[1]]@CIUpper.Est
      egg_est[i,] = m[[1]]@Estimate
      ivw_CI_lower[i,] = m[[2]]@CILower
      ivw_CI_upper[i,] = m[[2]]@CIUpper
      ivw_est[i,] = m[[2]]@Estimate
    }
    
    for (j in 1:3){
      egg_coverage_count[j] = egg_coverage_count[j] - sum(egg_CI_lower[,j] > beta[j]) - sum(egg_CI_upper[,j] < beta[j])
      ivw_coverage_count[j] = ivw_coverage_count[j] - sum(ivw_CI_lower[,j] > beta[j]) - sum(ivw_CI_upper[,j] < beta[j])
      egg_pwr_count[j] = sum(egg_CI_lower[,j] > 0) + sum(egg_CI_upper[,j] < 0)
      ivw_pwr_count[j] = sum(ivw_CI_lower[,j] > 0) + sum(ivw_CI_upper[,j] < 0)
    }
    
    results = rbind(results, data.frame(toString(alphamax),toString(measerror),
                                        mean(egg_est[,1]),sd(egg_est[,1]),egg_coverage_count[1]/trials,egg_pwr_count[1]/trials,
                                        mean(ivw_est[,1]),sd(ivw_est[,1]),ivw_coverage_count[1]/trials,ivw_pwr_count[1]/trials,
                                        mean(egg_est[,2]),sd(egg_est[,2]),egg_coverage_count[2]/trials,egg_pwr_count[2]/trials,
                                        mean(ivw_est[,2]),sd(ivw_est[,2]),ivw_coverage_count[2]/trials,ivw_pwr_count[2]/trials,
                                        mean(egg_est[,3]),sd(egg_est[,3]),egg_coverage_count[3]/trials,egg_pwr_count[3]/trials,
                                        mean(ivw_est[,3]),sd(ivw_est[,3]),ivw_coverage_count[3]/trials,ivw_pwr_count[3]/trials))
    
  }
}
