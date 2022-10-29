#' Compare simulated with real data
#' 
#' Violin plot to compare the mean and variance of simulated and real data.
#' @param sim_data simulated data
#' @param obs_datas observed data
#' @export
violin_plot<-function(sim_data, obs_data){
  mean_sim <- colMeans(sim_data)
  mean_sim_data<-data.frame("dataset" = rep("sim", ncol(sim_data)), "measure" = rep("mean", ncol(sim_data)), "value" = mean_sim)
  
  var_sim <- apply(sim_data,2,var)
  var_sim_data<-data.frame("dataset" = rep("sim", ncol(sim_data)), "measure" = rep("var", ncol(sim_data)), "value" = var_sim)
  
  mean_obs <- colMeans(obs_data)
  mean_obs_data<-data.frame("dataset" = rep("obs", ncol(obs_data)), "measure" = rep("mean", ncol(obs_data)), "value" = mean_obs)
  
  var_obs <- unname(apply(obs_data,2,var))
  var_obs_data <- data.frame("dataset" = rep("obs", ncol(obs_data)), "measure" = rep("var", ncol(obs_data)), "value" = mean_obs)
  
  stat<-bind_rows(mean_sim_data,var_sim_data, mean_obs_data,var_obs_data)
  
  ggplot(stat, aes(x = dataset, y = value)) +
    scale_y_log10()+
    geom_violin(scale = 'width', trim = TRUE) +
    facet_wrap(~measure, scales = "free", ncol = 3)
}