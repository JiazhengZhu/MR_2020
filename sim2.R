library(MendelianRandomization)
library(progress)
library(purrr)
library(doParallel)
library(foreach)

registerDoParallel(detectCores())

get_population_data = function(n_samples, n_snp, prob_snp, alpha_lim, beta_vec, xerror_vec, yerror, measerror_vec) {
  alpha_mat = cbind(runif(n_snp, min=alpha_lim[1], max = alpha_lim[2]), runif(n_snp, min=alpha_lim[1], max = alpha_lim[2]))
  n_exposures = length(beta_vec)
  g = matrix(rbinom(n_snp*n_samples,2,prob_snp), ncol = n_snp, nrow = n_samples)
  exposures = g%*%alpha_mat + t(t(matrix(rnorm(n_samples*n_exposures), ncol = n_exposures, nrow = n_samples))*xerror_vec)
  
  {####
    exposures[,2] = exposures[,1] + exposures[,2]
  }####
  
  y = exposures%*%beta_vec + rnorm(n_samples)*yerror
  
  exposures = exposures + t(t(matrix(rnorm(n_samples*n_exposures), ncol = n_exposures, nrow = n_samples))*measerror_vec)
  
  return(list(g,exposures,y))
}

get_firststage_values = function(l, choice = "all") {
  g = l[[1]]
  exposures = l[[2]]
  n_exposures = ncol(exposures)
  y = l[[3]]
  #G-Y association
  if (choice == "all" | choice == "gy") {
    by = rep(0, n_snp)
    byse = rep(0, n_snp)
    for (i in 1:n_snp){
      fit = lm(y~g[,i])
      by[i] = coef(fit)[[2]]
      byse[i] = coef(summary(fit))[, "Std. Error"][[2]]
    }
  }
  
  #G-X association
  bx = matrix(0, nrow = n_snp, ncol = n_exposures)
  bxse = matrix(0, nrow = n_snp, ncol = n_exposures)
  if (choice == "all" | choice == "gx") {
    for (j in 1:n_exposures){
      for (i in 1:n_snp){
        fit = lm(exposures[,j]~g[,i])
        bx[i,j] = coef(fit)[[2]]
        bxse[i,j] = coef(summary(fit))[, "Std. Error"][[2]]
      }
    }
  }
  
  if (choice == "all"){
    return(list(bx,bxse,by,byse))
  } else if (choice == "gx") {
    return(list(bx,bxse))
  } else if (choice == "gy") {
    return(list(by,byse))
  }
}

get_mv_models = function(l){
  bx = l[[1]]
  bxse = l[[2]]
  by = l[[3]]
  byse = l[[4]]
  mrmv = mr_mvinput(bx=bx,
                    bxse=bxse,
                    by=as.double(by),
                    byse=as.double(byse))
  return(list(mr_mvegger(mrmv), mr_mvivw(mrmv)))
}

sim_repeats = function(repeats, n_samples, n_snp, prob_snp, alpha_lim, beta_vec, xerror_vec, yerror, measerror_vec){
  n_exposures = length(beta_vec)
  egg_CI_lower = matrix(0, nrow = repeats, ncol = n_exposures)
  egg_CI_upper = matrix(0, nrow = repeats, ncol = n_exposures)
  egg_est = matrix(0, nrow = repeats, ncol = n_exposures)
  egg_int_CI_lower = matrix(0, nrow = repeats, ncol = 1)
  egg_int_CI_upper = matrix(0, nrow = repeats, ncol = 1)
  egg_int_est = matrix(0, nrow = repeats, ncol = 1)
  ivw_CI_lower = matrix(0, nrow = repeats, ncol = n_exposures)
  ivw_CI_upper = matrix(0, nrow = repeats, ncol = n_exposures)
  ivw_est = matrix(0, nrow = repeats, ncol = n_exposures)
  
  egg_coverage_count = ivw_coverage_count = rep(repeats,n_exposures)
  egg_pwr_count = ivw_pwr_count = rep(repeats,n_exposures)
  egg_int_coverage_count = 0
  
  pb <- progress_bar$new(format = "[:bar] :percent eta: :eta",
                         total = repeats)
  
  
  for (i in 1:repeats){
    pb$tick()
    #G-Y
    dat = get_population_data(n_samples, n_snp, prob_snp, alpha_lim, beta_vec, xerror_vec, yerror, measerror_vec)
    res = get_firststage_values(dat, choice = "gy")
    by = res[[1]]
    byse = res[[2]]
    #G-X
    dat = get_population_data(n_samples, n_snp, prob_snp, alpha_lim, beta_vec, xerror_vec, yerror, measerror_vec)
    res = get_firststage_values(dat, choice = "gx")
    bx = res[[1]]
    bxse = res[[2]]
    
    # dat = get_population_data(n_samples, n_snp, prob_snp, alpha_lim, beta_vec, xerror_vec, yerror, measerror_vec)
    # res = get_firststage_values(dat, choice = "gx")
    # bx[,2] = res[[1]][,2]
    # bxse[,2] = res[[2]][,2]
    # 
    
    #fit
    m = get_mv_models(list(bx,bxse,by,byse))
    
    egg_CI_lower[i,] = m[[1]]@CILower.Est
    egg_CI_upper[i,] = m[[1]]@CIUpper.Est
    egg_est[i,] = m[[1]]@Estimate
    egg_int_CI_lower[i,] = m[[1]]@CILower.Int
    egg_int_CI_upper[i,] = m[[1]]@CIUpper.Int
    egg_int_est[i,] = m[[1]]@Intercept
    ivw_CI_lower[i,] = m[[2]]@CILower
    ivw_CI_upper[i,] = m[[2]]@CIUpper
    ivw_est[i,] = m[[2]]@Estimate
  }
  for (j in 1:n_exposures){
    egg_coverage_count[j] = egg_coverage_count[j] - sum(egg_CI_lower[,j] > beta_vec[j]) - sum(egg_CI_upper[,j] < beta_vec[j])
    ivw_coverage_count[j] = ivw_coverage_count[j] - sum(ivw_CI_lower[,j] > beta_vec[j]) - sum(ivw_CI_upper[,j] < beta_vec[j])
    egg_pwr_count[j] = sum(egg_CI_lower[,j] > 0) + sum(egg_CI_upper[,j] < 0)
    ivw_pwr_count[j] = sum(ivw_CI_lower[,j] > 0) + sum(ivw_CI_upper[,j] < 0)
  }
  egg_int_coverage_count = sum(egg_int_CI_lower[,1] > 0) + sum(egg_int_CI_upper[,1] < 0)
  result = data.frame(toString(alpha_lim), toString(beta_vec), toString(measerror_vec))
  for (k in 1:n_exposures){
    result = cbind(result, data.frame(mean(egg_est[,k]),sd(egg_est[,k]),egg_coverage_count[k]/repeats,egg_pwr_count[k]/repeats,
                                      mean(ivw_est[,k]),sd(ivw_est[,k]),ivw_coverage_count[k]/repeats,ivw_pwr_count[k]/repeats))
  }
  result = cbind(result, data.frame(mean(egg_int_est),egg_int_coverage_count/repeats))
  return(result)
}

load_sim_params = function(rp, n_sam, nn_snp, p_snp, xe_vec, ye) {
  wrapped = function(ll){
    sim_repeats(repeats = rp, n_samples = n_sam, n_snp = nn_snp, prob_snp = p_snp, alpha_lim = ll[[1]], beta_vec = ll[[2]], xerror_vec = xe_vec, yerror = ye, measerror_vec = ll[[3]])
  }
  return(wrapped)
}

repeats = 2000
n_samples = 2000
n_snp = 20
prob_snp = 0.5
alpha_lim = list(c(0.6,0.8), c(0.4,0.6), c(0.2,0.4))
beta_vec = list(c(2,0.0),c(1,0.0), c(2,1.0), c(2, -1))
xerror_vec = c(1,1)
yerror = 1
measerror_vec = list(c(0,0.0), c(1,0.0), c(2,0.0), c(3,0), c(1,1), c(2,2), c(3,3))

params_superlist = cross(list(alpha_lim, beta_vec, measerror_vec))

res_tab = foreach(i = params_superlist, .combine = rbind, .packages = c("MendelianRandomization", "progress")) %dopar% {
  f = load_sim_params(repeats, n_samples, n_snp, prob_snp, xerror_vec, yerror)
  return(f(i))
}

write.csv(res_tab, paste0("par_res", format(Sys.time(), "%d-%b-%Y %H.%M"), ".csv"))

# tab_dep = sim_repeats(repeats, n_samples, n_snp, prob_snp, alpha_lim, beta_vec, xerror_vec, yerror, measerror_vec)

stopImplicitCluster()

