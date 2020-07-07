library(plotly)
library(plyr)
library(MendelianRandomization)
library(progress)
library(doParallel)
library(foreach)
library(purrr)


cluster <- makeCluster(detectCores()-1)
registerDoParallel(cluster)
simulate = function(para_list){
  nsamples = para_list[[1]]
  alpha_max = para_list[[2]]
  beta = para_list[[3]]
  xerror = para_list[[4]]
  measerror = para_list[[5]]
  yerror = para_list[[6]]
  g = matrix(rbinom(30*nsamples,2,0.2), ncol = 30, nrow = nsamples)
  beta1 = beta[1]
  beta2 = beta[2]
  beta3 = beta[3]
  
  alpha1 = runif(30,min=0,max=alpha_max[1])
  alpha2 = runif(30,min=0,max=alpha_max[2])
  alpha3 = runif(30,min=0,max=alpha_max[3])
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
    by[i] = coefficients(fit)[1]
    byse[i] = coef(summary(fit))[, "Std. Error"][1]
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
      bx[j,i] = coefficients(fit)[1]
      bxse[j,i] = coef(summary(fit))[, "Std. Error"][1]
    }
  }

  #plug in the model
  mrmv = mr_mvinput(bx=bx,
                    bxse=bxse,
                    by=as.double(by),
                    byse=as.double(byse))
  
  m = list(mr_mvegger(mrmv), mr_mvivw(mrmv))
  
  egg_CI_lower = m[[1]]@CILower.Est
  egg_CI_upper = m[[1]]@CIUpper.Est
  egg_est = m[[1]]@Estimate
  egg_int_CI_lower = m[[1]]@CILower.Int
  egg_int_CI_upper = m[[1]]@CIUpper.Int
  egg_int_est = m[[1]]@Intercept
  ivw_CI_lower = m[[2]]@CILower
  ivw_CI_upper = m[[2]]@CIUpper
  ivw_est = m[[2]]@Estimate
  
  ans = list(beta, alpha_max, measerror,egg_est,egg_CI_lower,egg_CI_upper,egg_int_est,egg_int_CI_lower,egg_int_CI_upper,ivw_est,ivw_CI_lower,ivw_CI_upper)
  names(ans) <- c("beta","alpha_max","measerror","egg_est","egg_CI_lower","egg_CI_upper","egg_int_est","egg_int_CI_lower","egg_int_CI_upper","ivw_est","ivw_CI_lower","ivw_CI_upper")
  
  return(data.frame(t(ans)))
  
}

#true values
nsamples = list(10)
alpha_list = list(c(1,1,1),
                  c(0.8,0.8,0.8),
                  c(0.4,0.4,0.4),
                  c(0.2,0.2,0.2))
beta = list(c(0.2,-0.6,0.0))
xerror = list(c(0.2,0.2,0.2))
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
yerror = list(0.1)

params_superlist = cross(list(nsamples, alpha_list, beta, xerror, measerror_list, yerror))

trials = 20

generate_data = function(trials, para_list){
  t = foreach(i=1:trials, .combine = rbind) %dopar% {
    simulate(para_list)
  }
  return(t)
}

collect_results = function(t){
  trials = nrow(t)
  nb = length(t$beta[[1]])
  egg_est = matrix(0, nrow = trials, ncol = nb)
  ivw_est = matrix(0, nrow = trials, ncol = nb)
  egg_int_est = matrix(0, nrow = trials, ncol = 1)
  egg_coverage_count = ivw_coverage_count = rep(trials,nb)
  egg_pwr_count = ivw_pwr_count = rep(trials,nb)
  egg_int_coverage_count = 0
  
  for (i in 1:trials){
    egg_est[i, ] = t$egg_est[[i]]
    ivw_est[i, ] = t$ivw_est[[i]]
    egg_int_est[i] = t$egg_int_est[[i]]
    for (j in 1:nb){
      if (t$egg_CI_lower[[i]][j] > t$beta[[i]][j] | t$egg_CI_upper[[i]][j] < t$beta[[i]][j]) {
      egg_coverage_count[j] = egg_coverage_count[j] -1
      }
      if (t$ivw_CI_lower[[i]][j] > t$beta[[i]][j] | t$ivw_CI_upper[[i]][j] < t$beta[[i]][j]) {
        ivw_coverage_count[j] = ivw_coverage_count[j] -1
      }
      if (t$egg_CI_lower[[i]][j] > 0 | t$egg_CI_upper[[i]][j] < 0) {
        egg_pwr_count[j] = egg_pwr_count[j] -1
      }
      if (t$ivw_CI_lower[[i]][j] > 0 | t$ivw_CI_upper[[i]][j] < 0) {
        ivw_pwr_count[j] = ivw_pwr_count[j] -1
      }
    }
    if (t$egg_int_CI_lower[[i]]>0 | t$egg_int_CI_lower[[i]]<0) {
      egg_int_coverage_count = egg_int_coverage_count + 1
    }
  }
  return(data.frame(toString(t$alpha_max),toString(t$measerror),
                                      mean(egg_est[,1]),sd(egg_est[,1]),egg_coverage_count[1]/trials,egg_pwr_count[1]/trials,
                                      mean(ivw_est[,1]),sd(ivw_est[,1]),ivw_coverage_count[1]/trials,ivw_pwr_count[1]/trials,
                                      mean(egg_est[,2]),sd(egg_est[,2]),egg_coverage_count[2]/trials,egg_pwr_count[2]/trials,
                                      mean(ivw_est[,2]),sd(ivw_est[,2]),ivw_coverage_count[2]/trials,ivw_pwr_count[2]/trials,
                                      mean(egg_est[,3]),sd(egg_est[,3]),egg_coverage_count[3]/trials,egg_pwr_count[3]/trials,
                                      mean(ivw_est[,3]),sd(ivw_est[,3]),ivw_coverage_count[3]/trials,ivw_pwr_count[3]/trials,
                                      mean(egg_int_est),sd(egg_int_est),egg_int_coverage_count/trials))
}

wrapped_sim = function(trials, para_list){
  t = generate_data(trials, para_list)
  return(collect_results(t))
}

results = data.frame()
for (i in params_superlist) {
  results = rbind(results, wrapped_sim(trials,i))
}

write.csv(results, file = "results2.csv")

stopCluster(cluster)

