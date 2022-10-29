library(tidyverse)
library(cmdstanr)
library(dplyr)
library(MASS)
library(MCMCpack)
library(markovchain)
library(tidybayes)

##fitting method for cluster
#' @export
cal_fit <- function(data_list) {
  ##npp_model <- cmdstan_model("./lmnsim/R/lmn_state_fitting.stan")
  npp_model <- cmdstan_model(system.file("lmn_state_fitting.stan", package = "lmnsim"))
  fit <-  npp_model$variational(data_list, iter = 2e4)
  a = fit$summary(variables = c("theta"), "mean")
  b = fit$summary(variables = c("sigma"), "mean")
  beta = fit$summary(variables = c("beta"),"mean")
  beta_matrix = matrix(nrow = data_list$n_clust, ncol = data_list$n_species)
  for (i in 1:(data_list$n_clust)) {
    for (j in 1:data_list$n_species) {
      beta_matrix[i, j] = beta$mean[beta$variable == paste("beta[", i, ",", j, "]", sep = "")]
    }
  }
  theta = data.frame()
  for (i in 1:(data_list$n_steps * data_list$n_person)) {
    for (j in 1:data_list$n_clust) {
      theta[i, j] = a$mean[a$variable == paste("theta[", i, ",", j, "]", sep = "")]
    }
  }
  sequence = data.frame()
  k = 1
  for (i in 1:data_list$n_person) {
    for (j in 1:data_list$n_steps) {
      sequence[i, j] =  which.max(theta[k, ])
      k = k + 1
    }
  }
  mcFit <- markovchainFit(data = sequence[1, ])
  fit_matrix = mcFit[["estimate"]]@transitionMatrix
  ##fix the problem of unvisited states
  if (nrow(fit_matrix) != data_list$n_clust || ncol(fit_matrix) != data_list$n_clust) {
    updated_matrix = matrix(nrow = data_list$n_clust, ncol = data_list$n_clust)
    for (j in 1:data_list$n_clust) {
      if (!(j %in% sequence)) {
        updated_matrix[1:ncol(updated_matrix), j] = 0
        updated_matrix[j, 1:ncol(updated_matrix)] = 1 / ncol(updated_matrix)
      }
    }
    for (j in rownames(fit_matrix)) {
      for (k in colnames(fit_matrix)) {
        updated_matrix[as.numeric(j), as.numeric(k)] = fit_matrix[j, k]
      }
    }
    return(list(sigma = b, matrix = updated_matrix, centers = beta_matrix))
  }
  return(list(sigma = b, matrix = fit_matrix, centers = beta_matrix))
}

#' Simulate microbiome sample data through cluster model
#' 
#' Simulate microbiome abundance count data with user-specified number of samples, species and sequencing depth
#' @param n_samples number of samples
#' @param n_species number of species
#' @param n_depth sequence depth
#' @param Sigma correlation matrix for muti-nomal distribution
#' @return the simulated data
#' @examples
#' sim_clust_obs(initial_state = 1, n_clust =4, n_person = 1, obs=preg_data)
#' @export
sim_clust_obs<-function(prob_matrix, initial_state, n_clust,n_person, sigma1, n_species, n_depth, obs){
  n_species = ncol(obs)-1
  n_steps = nrow(obs)/n_person
  n_depth = round(mean(rowSums(obs)),0)
  data_list <- list(n_species = n_species, n_steps = n_steps, n_clust = n_clust, y = obs, n_person = n_person)
  result = cal_fit(data_list)
  sigma1 = result$sigma
  sigma1 = sigma1$mean
  prob_matrix = result$matrix
  centers = result$centers
  states=markov_sample(prob_matrix, (n_steps*n_person), initial_state)
  data=matrix(nrow=(n_steps*n_person), ncol=n_species+1)
  for (i in 1:(n_steps*n_person)){
    x = centers[states[i],]
    data[i,]=rmultinom(n=1, size = n_depth, prob = phi_inverse(x))
  }
  list(states=states, data=data)
}

##fitting method for perturbation
model_fit<-function(data_list){
  model<-cmdstan_model("./lmnsim/R/lmn_perturbation_fitting.stan")
  fit <- model$variational(data_list, iter = 2e4)
  sigma = fit$summary(variables = c("sigma"), "mean")
  a = fit$summary(variables = c("beta"), "mean")
  beta = data.frame()
  for(i in 1:2){
    for (j in 1:data_list$n_species){
      beta[i,j] = a$mean[a$variable==paste("beta[",i,",", j, "]", sep = "")]
    } 
  }
  list(sigma = sigma$mean, beta = beta)
}

#' Simulate microbiome sample data through perturbation model
#' 
#' Simulate microbiome abundance count data with user-specified number of samples, species and sequencing depth
#' @param alpha the perturbation intensity, 0 means no perturbed, 1 means fully perturbed
#' @param obs observed microbiome data
#' @param n_depth sequence depth
#' @return the simulated data
#' @examples
#' sim_perturb_obs(alpha = c(0,0,0.8,0.8,0.8,0,0,0,0,0), obs = data, n_depth = 2000)
#' @export
sim_perturb_obs<-function(alpha, obs, n_depth = NULL){
  n_species<-ncol(obs)-1
  if(is.null(n_depth)){
    n_depth = round(mean(rowSums(obs)),0) 
  }
  data_list<-list(alpha=alpha, n_species=ncol(obs)-1, n_timepoints = length(alpha), y = obs)
  result = model_fit(data_list)
  beta = result$beta

  data = sim_ts_perturb_2(alpha = alpha, n_species = n_species, n_depth = n_depth, beta = beta)
  return (data)
}


