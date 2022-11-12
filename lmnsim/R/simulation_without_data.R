library(compositions)
library(tidyverse)
library(markovchain)
library(MCMCpack)
library(dplyr)
library(MASS)

phi_inverse <- function(mu) {
  sum_exp <- sum(exp(mu))
  c(exp(mu) / (1 + sum_exp), 1 / (1 + sum_exp))
}

#' Simulate microbiome sample data
#'
#' Simulate microbiome abundance count data with user-specified number of samples, species and sequencing depth
#' @param n_samples number of samples
#' @param n_species number of species
#' @param n_depth sequence depth
#' @param Sigma correlation matrix for muti-nomal distribution
#' @return the simulated data
#' @examples
#' simulate_samples(100,100,1000)
#' @export
sim_samples <- function(n_samples, n_species, n_depth, Sigma = NULL) {
  if (is.null(Sigma)) {
    Sigma <- matrix(0.1, n_species - 1, n_species - 1)
    diag(Sigma) <- 1
  }

  result <- list()
  x <- MASS::mvrnorm(1, rep(0, n_species - 1), Sigma)
  p <- phi_inverse(x)

  for (i in seq_len(n_samples)) {
    result[[i]] <- rmultinom(1, n_depth, p) %>%
      t() %>%
      as.data.frame()
  }

  bind_rows(result)
}

#' Simulate microbiome time series data
#'
#' simulate microbiome longitudinal count data for one sample
#' @param endtime the number of time points of the sample
#' @param n_species number of species
#' @param n_depth sequence depth
#' @param Sigma correlation matrix for muti-nomal distribution
#' @return the simulated data
#' @examples
#' simulate_ts(20,100,1000)
#' @export
sim_ts <- function(endtime, n_species, n_depth, Sigma = NULL) {
  if (is.null(Sigma)) {
    Sigma <- matrix(0.1, n_species - 1, n_species - 1)
    diag(Sigma) <- 1
  }

  result <- list()
  x <- MASS::mvrnorm(1, rep(0, n_species - 1), Sigma)
  for (i in seq_len(endtime)) {
    e <- MASS::mvrnorm(1, rep(0, n_species - 1), 0.001 * Sigma)
    x <- 0.999*x + e

    result[[i]] <- rmultinom(1, n_depth, phi_inverse(x)) %>%
      t() %>%
      as.data.frame()
  }

  bind_rows(result)
}

#' Simulate microbiome perturbation data
#'
#' To simulate microbiome longitudinal count data for perturbation experiments with user defined starting and ending time point, effect size and duration of perturbation
#' @param endtime the number of time points of the sample
#' @param n_species number of species
#' @param n_depth sequence depth
#' @param Sigma correlation matrix for muti-nomal distribution
#' @param beta effect of perturbation
#' @param t1 time of reaching the maximum of perturbation
#' @param t2 time of ending the maximum of perturbation
#' @param time the start time for perturbation
#' @param species perturbed species
#' @param duration the duration of perturbation
#' @return the simulated data
#' @examples
#' sim_ts_perturb(endtime=30, n_species=30, 1000, Sigma=NULL, beta=1,species=c(1,2,3),time=4, t1=8,t2=12,duration = 16)
#' @export
sim_ts_perturb <- function(endtime, n_species, n_depth, Sigma = NULL, beta, t1,t2, time, species, duration) {
  if (is.null(Sigma)) {
    Sigma <- matrix(0, n_species - 1, n_species - 1)
    diag(Sigma) <- 1
  }

  result <- list()
  x <- MASS::mvrnorm(1, rep(0, n_species - 1), Sigma)
  for (i in seq_len(endtime)) {
    if ((i>=time)&&(i<=(time+duration))){
      mean<-rep(0, n_species-1)
      func<-function(t1, t2, beta){
        effect=0
        if (i<t1){
          effect= beta*(i-time)
        }
        if(i>=t1 && i<=t2){
          effect = beta*(t1-time)
        }
        if (i>t2){
          effect = beta*(t1-time)-beta*(i-t2)
        }
        return (effect)
      }
      mean[species]<-func(t1, t2, beta)
      e <- MASS::mvrnorm(1, mean, 0.001 * Sigma)
      x <- 0.5*x + e
    }
    else{
      e <- MASS::mvrnorm(1, rep(0, n_species - 1), 0.001 * Sigma)
      x <- 0.5*x + e
    }
    result[[i]] <- rmultinom(1, n_depth, phi_inverse(x)) %>%
      t() %>%
      as.data.frame()
  }

  bind_rows(result)
}

#' Simulate microbiome perturbation data using clusted model
#'
#' To simulate microbiome longitudinal count data for perturbation experiments with user defined alpha, effect, and centroids of perturbed/unperturbed data
#' @param alpha the perturbation intensity, 0 means no perturbed, 1 means fully perturbed
#' @param n_species number of species
#' @param Sigma variance of different speicies in the same cluster
#' @param n_depth sequence depth
#' @param beta effect of perturbation
#' @param effect the effect of perturbation
#' @param beta the centroids of perturbed/unperturbed cluster
#' @return the simulated data
#' @examples
#' stan_sim(alpha=alpha, n_depth = 1000, n_species = 30, Sigma=NULL, effect=0)
#' @export
sim_ts_perturb_2<-function(alpha, n_species = 30, Sigma=NULL, n_depth, effect, beta = NULL){
  if (is.null(Sigma)) {
    Sigma <- matrix(0.1, n_species - 1, n_species - 1)
    diag(Sigma) <- 1
  }

  if(is.null(beta)){
    beta = matrix(nrow = 2, ncol = n_species-1)
    beta[1, ] = mvrnorm(n=1,rep(0, n_species - 1), Sigma)
    beta[2, ] = mvrnorm(n=1,c(rep(effect,3),rep(0, n_species - 4)), Sigma)
  }

  result <- list()
  for (i in 1:length(alpha)){

    x = c(t(beta) %*%  as.matrix(c(alpha[i], 1-alpha[i])))
    result[[i]] = rmultinom(1, n_depth, phi_inverse(x))%>%
      t() %>%
      as.data.frame()
  }

  as.matrix(bind_rows(result))
}

#' simulate the markov chain from transition matrix
markov_sample<-function(prob_matrix, n_steps, initial_state){
  result<-vector(length=n_steps)

  result[1]=initial_state
  for (i in 1:n_steps){

    result[i+1]=sample(x=1:nrow(prob_matrix),size=1, prob=prob_matrix[result[i],])
  }

  return(result)
}


#' Simulate microbiome clustered data
#'
#' To simulate microbiome longitudinal count data with clustering and markov chain based on user input transition matrix
#' @param prob_matrix transition matrix of markov chain
#' @param initial_state the initial state of the markov chain
#' @param n_steps number of time points
#' @param Sigma1 correlation matrix for centers
#' @param Sigma2 correlation matrix for the variation around centers
#' @param n_species number of species
#' @param n_depth sequence depth
#' @return the simulated data
#' @examples
#' sim_ts_cluster(prob_matrix=matrix(c(0.70, 0.2,0.1,0.3,0.4,0.3,0,0.4,0.6), nrow=3, byrow=TRUE), n_steps=25, initial_state = 3, sigma1=0.1, sigma2=0.01, n_species=100, n_depth=1000)
#' @export
sim_ts_cluster<-function(prob_matrix, initial_state , n_steps, sigma1, sigma2, n_species, n_depth){
  centers=mvrnorm(n=nrow(prob_matrix), rep(0, n_species-1), sigma1*diag(n_species-1))
  states=markov_sample(prob_matrix, n_steps, initial_state)
  data=matrix(nrow=n_steps, ncol=n_species)
  for (i in 1:n_steps){
    e = mvrnorm(n=1, rep(0, n_species-1), sigma2*diag(n_species-1))
    x = centers[states[i],] + e
    data[i,]=rmultinom(n=1, size = n_depth, prob = phi_inverse(x))
  }
  list(states=states, data=data)
}
