library(rlist)
#' Alpha sequence generator
#' 
#' The function generate alpha (perturbation intensity sequence) which can be used in the simulator 
#' @param n_timepoints number of time points
#' @param interval the interval of perturbation (interval[1]-interval[2]): increase the permutation intensity.(interval[2]-interval[3]): achieve maximized perturbation intensity. (interval[3]-interval[4]): decrease the permutation intensity.
#' @param signal_size the maximum perturbation intensity (between 0 and 1)
#' @return  alpha the perturbation intensity seuqence, 0 means no perturbed, 1 means fully perturbed
#' @export
make_alphas <- function(n_timepoints, intervals, signal_size) {
  alpha <- rep(0, n_timepoints)
  alpha[intervals[1]:intervals[2]] <- seq(0, signal_size, length.out = intervals[2] - intervals[1] + 1)
  alpha[intervals[2]:intervals[3]] <- signal_size
  alpha[intervals[3]:intervals[4]] <- seq(signal_size, 0, length.out = intervals[4] - intervals[3] + 1)
  alpha
}

#' configurations generator 
#' 
#' The function generate all configurations used in the power test 
#' @param ns the number of samples
#' @param n_depth sequence depth
#' @param alpha the maximum perturbation intensity (between 0 and 1)
#' @param interval the interval of perturbation (interval[1]-interval[2]): increase the permutation intensity.(interval[2]-interval[3]): achieve maximized perturbation intensity. (interval[3]-interval[4]): decrease the permutation intensity.
#' @param n_timepoints number of time points
#' @return configurations the configurations used in the power test
#' @export 
make_configurations <- function(ns, n_depth, alpha, intervals = NULL, n_timepoints = NULL) {
  if (is.null(intervals)) {
    intervals <- c(5, 8, 14, 18)
  }
  if (is.null(n_timepoints)) {
    n_timepoints <- 20
  }
  
  
  configurations <- cross(list(n = ns, n_depth = n_depth, alpha = alpha))
  for (i in seq_along(configurations)) {
    configurations[[i]][["alpha"]] <- make_alphas(n_timepoints, intervals, configurations[[i]][["alpha"]]) %>%
      round(3)
  }
  configurations
}

#' @export
get_centroids <- function(data_list){
  result = model_fit(data_list)
  result$beta
}


#' @export
modify_centroids_fraction <- function(centroids, fraction_perturbed = NULL) {
  # estimated centroids from pilot # get_centroids(...)
  # sort effects from highest to lowest (reorder columns)
  # find threshold for fraction_perturbed
  # all columns after fraction_perturbed -> replace with midpoint of two centroids
  if(is.null(fraction_perturbed)){
    fraction_perturbed = 0.5
  }
  
  difference = sort(rank(abs(centroids[2,] - centroids[1,])))
  for(i in 1:(ncol(centroids)*(1-fraction_perturbed))){
    column = as.integer(substring(names(difference[i]),2))
    centroids[1,column] = (centroids[1,column] + centroids[2,column])/2
    centroids[2,column] = centroids[1,column] 
  }
  centroids
}

#' @export
power_matrix_transform = function(matrix){
  stat = as.data.frame(matrix)
  for(i in 1:ncol(stat)){
    colnames(stat)[i] = paste("rep" , i, sep = "")
  }
  for(i in 1:length(configurations)){
    stat$n[i] = configurations[[i]]$n
    stat$seq_depth[i] = configurations[[i]]$seq_depth
    stat$alpha[i] = max(configurations[[i]]$alpha) ##signal_size, need to be changed if the graph changed.
  }
  pivot_longer(stat, cols = starts_with("rep"),names_to = "rep",values_to = "power")
}

#' @export
estimate_stat <-
  function(config,
           centroids,
           n_reps,
           simulator,
           sig_level,
           sim_opts) {
    ##Take average of samples, then do t-test
    statistics = matrix(nrow = length(configurations), ncol = n_reps)
    for (i in 1:length(configurations)) {
      for (j in 1:n_reps) {
        samples = data.frame()
        for (k in 1:configurations[[i]]$n) {
          samples = rbind(samples, simulator(
            alpha = configurations[[i]]$alpha,
            n_depth = configurations[[i]]$n_depth,
            beta = centroids
          ))
        }
        detected = 0
        truth = 0
        for (k in sim_opts) {
          perturbed = vector()
          unperturbed = vector()
          for (l in 1:(length(configurations[[i]]$alpha)*configurations[[i]]$n)) {
            if (configurations[[i]]$alpha[((l-1)%%length(configurations[[i]]$alpha))+1] == 0) {
              unperturbed = list.append(unperturbed, samples[l, k])
            }
            else{
              perturbed = list.append(perturbed, samples[l, k])
            }
          }
          truth = truth + 1
          if (t.test(perturbed, unperturbed)$p.value < sig_level) {
            detected = detected + 1
          }
        }
        statistics[i, j] = detected / truth
      }
    }
    statistics
  }