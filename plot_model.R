rm(list = ls())

source(file = "utils/functions.R"); source(file = "utils/retention_fit_top_up.R"); source(file = "utils/run_model_functions.R");
source(file = "read_data.R")

#########################
##### top up values #####
#########################

list2env(readRDS(file = "parameters/top_up_fits.rds"), envir = globalenv())

######################
##### model runs #####
######################

# model outputs from the run_model.R script
#list2env(readRDS(file = "results/model_runs.rds"), envir = globalenv())

list2env(readRDS(file = "results/mcmc_vals.rds"), envir = globalenv())
list2env(readRDS(file = "results/model_run_combs.rds"), envir = globalenv())
pred_prev_T <- readRDS(file = "results/pred_prev_T_b_mcmc.rds")
pred_prev_B <- readRDS(file = "results/pred_prev_B_b_mcmc.rds")

##########################  
##### data wrangling #####
##########################
# formatting the results
for(i in 1:length(pred_prev_B)){
  pred_prev_B[[i]]$date <- seq(start_date, by = "day", length.out = length(pred_prev_B[[i]]$timestep))
}

# for(i in 1:length(pred_prev_B_tu)){
#   pred_prev_B_tu[[i]]$date <- seq(start_date, by = "day", length.out = length(pred_prev_B_tu[[i]]$timestep))
# }
# 
# for(i in 1:length(pred_prev_B_tu_none)){
#   pred_prev_B_tu_none[[i]]$date <- seq(start_date, by = "day", length.out = length(pred_prev_B_tu_none[[i]]$timestep))
# }

for(i in 1:length(pred_prev_T)){
  pred_prev_T[[i]]$date <- seq(start_date, by = "day", length.out = length(pred_prev_T[[1]]$timestep))
}

# for(i in 1:length(pred_prev_T_tu)){
#   pred_prev_T_tu[[i]]$date <- seq(start_date, by = "day", length.out = length(pred_prev_T_tu[[1]]$timestep))
# }
# 
# for(i in 1:length(pred_prev_T_tu_none)){
#   pred_prev_T_tu_none[[i]]$date <- seq(start_date, by = "day", length.out = length(pred_prev_T_tu_none[[1]]$timestep))
# }

# data frames with the sampled prevalence values for plotting
df_T_prev <- na.omit(df_T[,c("Net", "date", "Malaria_prevalence", "Malaria_prevalence_l", "Malaria_prevalence_u")])
df_B_prev <- na.omit(df_B[,c("Net", "date", "Malaria_prevalence", "Malaria_prevalence_l", "Malaria_prevalence_u", "Location")])

#####################
##### functions #####
#####################

# calculates the prevalence for each unique model run and calculates the mean for repeated simulations
extract_prev <- function(Location_in,
                         Net_in,
                         top_up_name,
                         Country_in,
                         r_model,
                         int_Net_in,
                         top_up_net,
                         trial_net_cov_in,
                         vals_pred, 
                         pred_prev,
                         n_mcmc,
                         quant
                         ){
  
  # for the models with the top up of nets there is no trial_net_cov_in column
  # unless the nets are topped up with the trial nets or another type of net there is no top_up_net column
  
    inds <- which(vals_pred$Location_in == Location_in &
                  vals_pred$Net_in == Net_in &
                  vals_pred$top_up_name == top_up_name &
                  vals_pred$Country_in == Country_in &
                  vals_pred$r_model == r_model &
                  vals_pred$int_Net_in == int_Net_in &
                  vals_pred$top_up_net == top_up_net &
                  is.na(vals_pred$trial_net_cov_in) == is.na(trial_net_cov_in))
    
  if(length(inds) != n_mcmc){return(NA)}else{
    
    prev_mat <- if(Country_in == "Tanzania"){
      sapply(inds, function(i, p){summary_pfpr_0.5_14(p[[i]])}, p = pred_prev)
      } else if(Country_in == "Benin"){
        sapply(inds, function(i, p){summary_pfpr_0_100_all(p[[i]])}, p = pred_prev)
    }
    
    # returns the mean value from the stochastic simulations
    return(apply(prev_mat, 1, quantile, probs = c(quant)))
  }
}

extract_prev_arms <- function(r_in){
  arm_sims_T <- data.frame("t" = rep(pred_prev_T[[1]]$timestep, 3),
                           "date" = rep(pred_prev_T[[1]]$date, 3),
                           "net" = c(rep("IG2", sim_length),
                                     rep("PBO", sim_length),
                                     rep("Pyrethroid_only", sim_length)),
                           
                           "prev" = c(extract_prev(Location_in = "Misungwi",
                                                   Net_in = "IG2", 
                                                   top_up_name = "top_up_IG2_T",
                                                   Country_in = "Tanzania",
                                                   r_model = r_in,
                                                   int_Net_in = "IG2",
                                                   top_up_net = "Pyrethroid_only",
                                                   trial_net_cov_in = NA,
                                                   vals_pred = T_vals_pred, 
                                                   pred_prev = pred_prev_T,
                                                   n_mcmc = n_mcmc,
                                                   quant = 0.5),
                                      
                                      extract_prev(Location_in = "Misungwi",
                                                   Net_in = "PBO", 
                                                   top_up_name = "top_up_PBO_T",
                                                   Country_in = "Tanzania",
                                                   r_model = r_in,
                                                   int_Net_in = "PBO",
                                                   top_up_net = "Pyrethroid_only",
                                                   trial_net_cov_in = NA,
                                                   vals_pred = T_vals_pred, 
                                                   pred_prev = pred_prev_T,
                                                   n_mcmc = n_mcmc,
                                                   quant = 0.5),
                                      
                                      extract_prev(Location_in = "Misungwi",
                                                   Net_in = "Pyrethroid_only", 
                                                   top_up_name = "top_up_p_only_T",
                                                   Country_in = "Tanzania",
                                                   r_model = r_in,
                                                   int_Net_in = "Pyrethroid_only",
                                                   top_up_net = "Pyrethroid_only",
                                                   trial_net_cov_in = NA,
                                                   vals_pred = T_vals_pred, 
                                                   pred_prev = pred_prev_T,
                                                   n_mcmc = n_mcmc,
                                                   quant = 0.5)
                                      ),
                           
                           "lower" = c(extract_prev(Location_in = "Misungwi",
                                                    Net_in = "IG2", 
                                                    top_up_name = "top_up_IG2_T",
                                                    Country_in = "Tanzania",
                                                    r_model = r_in,
                                                    int_Net_in = "IG2",
                                                    top_up_net = "Pyrethroid_only",
                                                    trial_net_cov_in = NA,
                                                    vals_pred = T_vals_pred, 
                                                    pred_prev = pred_prev_T,
                                                    n_mcmc = n_mcmc,
                                                    quant = 0.025),
                                       
                                       extract_prev(Location_in = "Misungwi",
                                                    Net_in = "PBO", 
                                                    top_up_name = "top_up_PBO_T",
                                                    Country_in = "Tanzania",
                                                    r_model = r_in,
                                                    int_Net_in = "PBO",
                                                    top_up_net = "Pyrethroid_only",
                                                    trial_net_cov_in = NA,
                                                    vals_pred = T_vals_pred, 
                                                    pred_prev = pred_prev_T,
                                                    n_mcmc = n_mcmc,
                                                    quant = 0.025),
                                       
                                       extract_prev(Location_in = "Misungwi",
                                                    Net_in = "Pyrethroid_only", 
                                                    top_up_name = "top_up_p_only_T",
                                                    Country_in = "Tanzania",
                                                    r_model = r_in,
                                                    int_Net_in = "Pyrethroid_only",
                                                    top_up_net = "Pyrethroid_only",
                                                    trial_net_cov_in = NA,
                                                    vals_pred = T_vals_pred, 
                                                    pred_prev = pred_prev_T,
                                                    n_mcmc = n_mcmc,
                                                    quant = 0.025)
                                       ),
                           
                           "upper" = c(extract_prev(Location_in = "Misungwi",
                                                    Net_in = "IG2", 
                                                    top_up_name = "top_up_IG2_T",
                                                    Country_in = "Tanzania",
                                                    r_model = r_in,
                                                    int_Net_in = "IG2",
                                                    top_up_net = "Pyrethroid_only",
                                                    trial_net_cov_in = NA,
                                                    vals_pred = T_vals_pred, 
                                                    pred_prev = pred_prev_T,
                                                    n_mcmc = n_mcmc,
                                                    quant = 0.975),
                                       
                                       extract_prev(Location_in = "Misungwi",
                                                    Net_in = "PBO", 
                                                    top_up_name = "top_up_PBO_T",
                                                    Country_in = "Tanzania",
                                                    r_model = r_in,
                                                    int_Net_in = "PBO",
                                                    top_up_net = "Pyrethroid_only",
                                                    trial_net_cov_in = NA,
                                                    vals_pred = T_vals_pred, 
                                                    pred_prev = pred_prev_T,
                                                    n_mcmc = n_mcmc,
                                                    quant = 0.975),
                                       
                                       extract_prev(Location_in = "Misungwi",
                                                    Net_in = "Pyrethroid_only", 
                                                    top_up_name = "top_up_p_only_T",
                                                    Country_in = "Tanzania",
                                                    r_model = r_in,
                                                    int_Net_in = "Pyrethroid_only",
                                                    top_up_net = "Pyrethroid_only",
                                                    trial_net_cov_in = NA,
                                                    vals_pred = T_vals_pred, 
                                                    pred_prev = pred_prev_T,
                                                    n_mcmc = n_mcmc,
                                                    quant = 0.975)
                                       )
                           )
  
  arm_sims_B <- data.frame("t" = rep(pred_prev_B[[1]]$timestep, 2),
                           "date" = rep(pred_prev_B[[1]]$date, 2),
                           "net" = c(rep("IG2", sim_length),
                                     rep("Pyrethroid_only", sim_length)),
                           "prev" = c(extract_prev(Location_in = "Zagnanado",
                                                     Net_in = "IG2", 
                                                     top_up_name = "top_up_IG2_B",
                                                     Country_in = "Benin",
                                                     r_model = r_in,
                                                     int_Net_in = "IG2",
                                                     top_up_net = "Pyrethroid_only",
                                                     trial_net_cov_in = NA,
                                                     vals_pred = B_vals_pred, 
                                                     pred_prev = pred_prev_B,
                                                     n_mcmc = n_mcmc,
                                                     quant = 0.5),
                                      extract_prev(Location_in = "Ouinhi",
                                                     Net_in = "Pyrethroid_only", 
                                                     top_up_name = "top_up_p_only_B",
                                                     Country_in = "Benin",
                                                     r_model = r_in,
                                                     int_Net_in = "Pyrethroid_only",
                                                     top_up_net = "Pyrethroid_only",
                                                     trial_net_cov_in = NA,
                                                     vals_pred = B_vals_pred, 
                                                     pred_prev = pred_prev_B,
                                                     n_mcmc = n_mcmc,
                                                     quant = 0.5)
                                      ),
                           "lower" = c(extract_prev(Location_in = "Zagnanado",
                                                    Net_in = "IG2", 
                                                    top_up_name = "top_up_IG2_B",
                                                    Country_in = "Benin",
                                                    r_model = r_in,
                                                    int_Net_in = "IG2",
                                                    top_up_net = "Pyrethroid_only",
                                                    trial_net_cov_in = NA,
                                                    vals_pred = B_vals_pred, 
                                                    pred_prev = pred_prev_B,
                                                    n_mcmc = n_mcmc,
                                                    quant = 0.025),
                                       extract_prev(Location_in = "Ouinhi",
                                                    Net_in = "Pyrethroid_only", 
                                                    top_up_name = "top_up_p_only_B",
                                                    Country_in = "Benin",
                                                    r_model = r_in,
                                                    int_Net_in = "Pyrethroid_only",
                                                    top_up_net = "Pyrethroid_only",
                                                    trial_net_cov_in = NA,
                                                    vals_pred = B_vals_pred, 
                                                    pred_prev = pred_prev_B,
                                                    n_mcmc = n_mcmc,
                                                    quant = 0.025)
                                       ),
                           "upper" = c(extract_prev(Location_in = "Zagnanado",
                                                    Net_in = "IG2", 
                                                    top_up_name = "top_up_IG2_B",
                                                    Country_in = "Benin",
                                                    r_model = r_in,
                                                    int_Net_in = "IG2",
                                                    top_up_net = "Pyrethroid_only",
                                                    trial_net_cov_in = NA,
                                                    vals_pred = B_vals_pred, 
                                                    pred_prev = pred_prev_B,
                                                    n_mcmc = n_mcmc,
                                                    quant = 0.975),
                                       extract_prev(Location_in = "Ouinhi",
                                                    Net_in = "Pyrethroid_only", 
                                                    top_up_name = "top_up_p_only_B",
                                                    Country_in = "Benin",
                                                    r_model = r_in,
                                                    int_Net_in = "Pyrethroid_only",
                                                    top_up_net = "Pyrethroid_only",
                                                    trial_net_cov_in = NA,
                                                    vals_pred = B_vals_pred, 
                                                    pred_prev = pred_prev_B,
                                                    n_mcmc = n_mcmc,
                                                    quant = 0.975)
                                       )
  )
  
  T_arm_plot <- ggplot() +
    geom_vline(xintercept = intervention_date_T, linetype = 2, linewidth = 1, alpha = 0.5) +
    geom_ribbon(data = arm_sims_T,
                aes(x = date, y = prev, ymin = upper, ymax = lower, 
                    group = net, fill = net), alpha = 0.25) +
    geom_line(data = arm_sims_T,
              aes(x = date, y = prev, 
                  col = net, group = net), linewidth = 1) +
    geom_pointrange(data = subset(df_T_prev, 
                                  Net == "PBO"), 
                    aes(x = date, y = Malaria_prevalence, ymin = Malaria_prevalence_l, ymax = Malaria_prevalence_u),
                    fill = "aquamarine", size = 0.9, shape = 21) +
    geom_pointrange(data = subset(df_T_prev, Net == "IG2"), 
                    aes(x = date, y = Malaria_prevalence, ymin = Malaria_prevalence_l, ymax = Malaria_prevalence_u),
                    fill = "darkgreen", size = 0.9, shape = 21) +
    geom_pointrange(data = subset(df_T_prev, Net == "Pyrethroid_only"), 
                    aes(x = date, y = Malaria_prevalence, ymin = Malaria_prevalence_l, ymax = Malaria_prevalence_u),
                    fill = "blue", size = 0.9, shape = 21) +
    scale_colour_manual(breaks = c("Pyrethroid_only", "PBO", "IG2"), values = c("blue", "aquamarine", "darkgreen"), name = "",
                        labels = c("pyrethroid-only", "pyrethroid-PBO", "pyrethroid-pyrrole")) +
    scale_fill_manual(breaks = c("Pyrethroid_only", "PBO", "IG2"), values = c("blue","aquamarine", "darkgreen"), name = "",
                      labels = c("pyrethroid-only", "pyrethroid-PBO", "pyrethroid-pyrrole")) +
    xlab("Year") +
    ylab("Malaria prevalence in children\naged 0.5 to 14 years old") +
    coord_cartesian(xlim = c(baseline_start_date_T - 365/2, intervention_date_T + 365*3),
                    ylim = c(0, 0.7)) +
    theme_classic() +
    ggtitle("Tanzania trial arms") +
    theme(text = element_text(size = 20),
          axis.text = element_text(size = 17),
          legend.text = element_text(size = 14),
          legend.title = element_text(size=14),
          legend.background = element_blank(),
          legend.position = c(0.85, 0.95)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1), breaks = seq(0, 0.7, 0.1))
  
  B_arm_plot <- ggplot() +
    geom_vline(xintercept = intervention_date_B, linetype = 2, linewidth = 1, alpha = 0.5) +
    geom_ribbon(data = arm_sims_B, 
                aes(x = date, ymin = upper, ymax = lower, fill = net, group = net), 
                alpha = 0.25) +
    geom_line(data = arm_sims_B,
              aes(x = date, y = prev, col = net), linewidth = 1) +
    # geom_point(data = subset(df_B, Location == "Cove"), 
    #            aes(x = date, y = Malaria_prevalence),
    #            fill = "skyblue", size = 3.5, shape = 21) +
    geom_pointrange(data = subset(df_B_prev, Location == "Zagnanado"), 
                    aes(x = date, y = Malaria_prevalence, ymin = Malaria_prevalence_l, ymax = Malaria_prevalence_u),
                    fill = "darkgreen", size = 0.9, shape = 22) +
    geom_pointrange(data = subset(df_B_prev, Location == "Ouinhi"), 
                    aes(x = date, y = Malaria_prevalence, ymin = Malaria_prevalence_l, ymax = Malaria_prevalence_u),
                    fill = "blue", size = 0.9, shape = 22) +
    scale_colour_manual(values = c("blue", "darkgreen"), name = "", breaks = c("Pyrethroid_only", "IG2"),
                        labels = c("pyrethroid-only", "pyrethroid-pyrrole")) +
    scale_fill_manual(values = c("blue", "darkgreen"), name = "", breaks = c("Pyrethroid_only", "IG2"),
                      labels = c("pyrethroid-only", "pyrethroid-pyrrole")) +
    xlab("Year") +
    ylab("Malaria prevalence in\npeople of all ages") +
    coord_cartesian(xlim = c(baseline_start_date_B - 365/2, intervention_date_B + 365*3),
                    ylim = c(0, 0.7)) +
    theme_classic() +
    ggtitle("Benin trial arms") +
    theme(text = element_text(size = 20),
          axis.text = element_text(size = 17),
          legend.text = element_text(size = 14),
          legend.title = element_text(size=14),
          legend.background = element_blank(),
          legend.position = c(0.85, 0.95)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1), breaks = seq(0, 0.7, 0.1))
  
  return(list("arm_sims_T" = arm_sims_T, "arm_sims_B" = arm_sims_B,
              "T_arm_plot" = T_arm_plot, "B_arm_plot" = B_arm_plot))
}

# pyrethroid-pyrrole arms with pyrethroid-pyrrole top up nets

extract_prev_diff_top_up <- function(r_in){
  T_tu_df <- data.frame("t" = rep(pred_prev_T[[1]]$timestep, 3),
                      "date" = rep(pred_prev_T[[1]]$date, 3),
                      "net" = c(rep("IG2", sim_length),
                                rep("IG2_tu", sim_length),
                                rep("None", sim_length)),
                      
                      "prev" = c(extract_prev(Location_in = "Misungwi",
                                              Net_in = "IG2", 
                                              top_up_name = "top_up_IG2_T",
                                              Country_in = "Tanzania",
                                              r_model = r_in,
                                              int_Net_in = "IG2",
                                              top_up_net = "Pyrethroid_only",
                                              trial_net_cov_in = NA,
                                              vals_pred = T_vals_pred, 
                                              pred_prev = pred_prev_T,
                                              n_mcmc = n_mcmc,
                                              quant = 0.5),
                                 extract_prev(Location_in = "Misungwi",
                                              Net_in = "IG2", 
                                              top_up_name = "top_up_IG2_T",
                                              Country_in = "Tanzania",
                                              r_model = r_in,
                                              int_Net_in = "IG2",
                                              top_up_net = "IG2",
                                              trial_net_cov_in = NA,
                                              vals_pred = T_vals_pred, 
                                              pred_prev = pred_prev_T,
                                              n_mcmc = n_mcmc,
                                              quant = 0.5),
                                 extract_prev(Location_in = "Misungwi",
                                              Net_in = "IG2", 
                                              top_up_name = "top_up_IG2_T",
                                              Country_in = "Tanzania",
                                              r_model = r_in,
                                              int_Net_in = "IG2",
                                              top_up_net = "none",
                                              trial_net_cov_in = 0,
                                              vals_pred = T_vals_pred, 
                                              pred_prev = pred_prev_T,
                                              n_mcmc = n_mcmc,
                                              quant = 0.5) # no top up
                      ),
                      
                      "lower" = c(extract_prev(Location_in = "Misungwi",
                                               Net_in = "IG2", 
                                               top_up_name = "top_up_IG2_T",
                                               Country_in = "Tanzania",
                                               r_model = r_in,
                                               int_Net_in = "IG2",
                                               top_up_net = "Pyrethroid_only",
                                               trial_net_cov_in = NA,
                                               vals_pred = T_vals_pred, 
                                               pred_prev = pred_prev_T,
                                               n_mcmc = n_mcmc,
                                               quant = 0.025),
                                  extract_prev(Location_in = "Misungwi",
                                               Net_in = "IG2", 
                                               top_up_name = "top_up_IG2_T",
                                               Country_in = "Tanzania",
                                               r_model = r_in,
                                               int_Net_in = "IG2",
                                               top_up_net = "IG2",
                                               trial_net_cov_in = NA,
                                               vals_pred = T_vals_pred, 
                                               pred_prev = pred_prev_T,
                                               n_mcmc = n_mcmc,
                                               quant = 0.025),
                                  extract_prev(Location_in = "Misungwi",
                                               Net_in = "IG2", 
                                               top_up_name = "top_up_IG2_T",
                                               Country_in = "Tanzania",
                                               r_model = r_in,
                                               int_Net_in = "IG2",
                                               top_up_net = "none",
                                               trial_net_cov_in = 0,
                                               vals_pred = T_vals_pred, 
                                               pred_prev = pred_prev_T,
                                               n_mcmc = n_mcmc,
                                               quant = 0.025)
                      ),
                      
                      "upper" = c(extract_prev(Location_in = "Misungwi",
                                               Net_in = "IG2", 
                                               top_up_name = "top_up_IG2_T",
                                               Country_in = "Tanzania",
                                               r_model = r_in,
                                               int_Net_in = "IG2",
                                               top_up_net = "Pyrethroid_only",
                                               trial_net_cov_in = NA,
                                               vals_pred = T_vals_pred, 
                                               pred_prev = pred_prev_T,
                                               n_mcmc = n_mcmc,
                                               quant = 0.975),
                                  extract_prev(Location_in = "Misungwi",
                                               Net_in = "IG2", 
                                               top_up_name = "top_up_IG2_T",
                                               Country_in = "Tanzania",
                                               r_model = r_in,
                                               int_Net_in = "IG2",
                                               top_up_net = "IG2",
                                               trial_net_cov_in = NA,
                                               vals_pred = T_vals_pred, 
                                               pred_prev = pred_prev_T,
                                               n_mcmc = n_mcmc,
                                               quant = 0.975),
                                  extract_prev(Location_in = "Misungwi",
                                               Net_in = "IG2", 
                                               top_up_name = "top_up_IG2_T",
                                               Country_in = "Tanzania",
                                               r_model = r_in,
                                               int_Net_in = "IG2",
                                               top_up_net = "none",
                                               trial_net_cov_in = 0,
                                               vals_pred = T_vals_pred, 
                                               pred_prev = pred_prev_T,
                                               n_mcmc = n_mcmc,
                                               quant = 0.975)
                      )
                      )
  
  B_tu_df <- data.frame("t" = rep(pred_prev_B[[1]]$timestep, 3),
                      "date" = rep(pred_prev_B[[1]]$date, 3),
                      "net" = c(rep("IG2", sim_length),
                                rep("IG2_tu", sim_length),
                                rep("None", sim_length)),
                      
                      "prev" = c(extract_prev(Location_in = "Zagnanado",
                                              Net_in = "IG2", 
                                              top_up_name = "top_up_IG2_B",
                                              Country_in = "Benin",
                                              r_model = r_in,
                                              int_Net_in = "IG2",
                                              top_up_net = "Pyrethroid_only",
                                              trial_net_cov_in = NA,
                                              vals_pred = B_vals_pred, 
                                              pred_prev = pred_prev_B,
                                              n_mcmc = n_mcmc,
                                              quant = 0.5),
                                 
                                 extract_prev(Location_in = "Zagnanado",
                                              Net_in = "IG2", 
                                              top_up_name = "top_up_IG2_B",
                                              Country_in = "Benin",
                                              r_model = r_in,
                                              int_Net_in = "IG2",
                                              top_up_net = "IG2",
                                              trial_net_cov_in = NA,
                                              vals_pred = B_vals_pred, 
                                              pred_prev = pred_prev_B,
                                              n_mcmc = n_mcmc,
                                              quant = 0.5),
                                 
                                 extract_prev(Location_in = "Zagnanado",
                                              Net_in = "IG2", 
                                              top_up_name = "top_up_IG2_B",
                                              Country_in = "Benin",
                                              r_model = r_in,
                                              int_Net_in = "IG2",
                                              top_up_net = "none",
                                              trial_net_cov_in = 0,
                                              vals_pred = B_vals_pred, 
                                              pred_prev = pred_prev_B,
                                              n_mcmc = n_mcmc,
                                              quant = 0.5)
                      ),
                      
                      "lower" = c(extract_prev(Location_in = "Zagnanado",
                                               Net_in = "IG2", 
                                               top_up_name = "top_up_IG2_B",
                                               Country_in = "Benin",
                                               r_model = r_in,
                                               int_Net_in = "IG2",
                                               top_up_net = "Pyrethroid_only",
                                               trial_net_cov_in = NA,
                                               vals_pred = B_vals_pred, 
                                               pred_prev = pred_prev_B,
                                               n_mcmc = n_mcmc,
                                               quant = 0.025),
                                  
                                  extract_prev(Location_in = "Zagnanado",
                                               Net_in = "IG2", 
                                               top_up_name = "top_up_IG2_B",
                                               Country_in = "Benin",
                                               r_model = r_in,
                                               int_Net_in = "IG2",
                                               top_up_net = "IG2",
                                               trial_net_cov_in = NA,
                                               vals_pred = B_vals_pred, 
                                               pred_prev = pred_prev_B,
                                               n_mcmc = n_mcmc,
                                               quant = 0.025),
                                  
                                  extract_prev(Location_in = "Zagnanado",
                                               Net_in = "IG2", 
                                               top_up_name = "top_up_IG2_B",
                                               Country_in = "Benin",
                                               r_model = r_in,
                                               int_Net_in = "IG2",
                                               top_up_net = "none",
                                               trial_net_cov_in = 0,
                                               vals_pred = B_vals_pred, 
                                               pred_prev = pred_prev_B,
                                               n_mcmc = n_mcmc,
                                               quant = 0.025)
                      ),
                      
                      "upper" = c(extract_prev(Location_in = "Zagnanado",
                                               Net_in = "IG2", 
                                               top_up_name = "top_up_IG2_B",
                                               Country_in = "Benin",
                                               r_model = r_in,
                                               int_Net_in = "IG2",
                                               top_up_net = "Pyrethroid_only",
                                               trial_net_cov_in = NA,
                                               vals_pred = B_vals_pred, 
                                               pred_prev = pred_prev_B,
                                               n_mcmc = n_mcmc,
                                               quant = 0.975),
                                  
                                  extract_prev(Location_in = "Zagnanado",
                                               Net_in = "IG2", 
                                               top_up_name = "top_up_IG2_B",
                                               Country_in = "Benin",
                                               r_model = r_in,
                                               int_Net_in = "IG2",
                                               top_up_net = "IG2",
                                               trial_net_cov_in = NA,
                                               vals_pred = B_vals_pred, 
                                               pred_prev = pred_prev_B,
                                               n_mcmc = n_mcmc,
                                               quant = 0.975),
                                  
                                  extract_prev(Location_in = "Zagnanado",
                                               Net_in = "IG2", 
                                               top_up_name = "top_up_IG2_B",
                                               Country_in = "Benin",
                                               r_model = r_in,
                                               int_Net_in = "IG2",
                                               top_up_net = "none",
                                               trial_net_cov_in = 0,
                                               vals_pred = B_vals_pred, 
                                               pred_prev = pred_prev_B,
                                               n_mcmc = n_mcmc,
                                               quant = 0.975)
                      )
)
  
  cf_plot_T <- ggplot() +
  geom_ribbon(data = subset(T_tu_df, date <= intervention_date_T),
              aes(x = date, ymin = upper, ymax = lower, fill = net, group = net),
              alpha = 0.1) +
  
  geom_line(data = subset(T_tu_df, date <= intervention_date_T),
            aes(x = date, y = prev, col = net, group = net),
            alpha = 0.1, linewidth = 0.2) +
  
  geom_ribbon(data = subset(T_tu_df, date >= intervention_date_T),
              aes(x = date, ymin = upper, ymax = lower, fill = net, group = net),
              alpha = 0.25) +
  
  geom_line(data = subset(T_tu_df, date >= intervention_date_T),
            aes(x = date, y = prev, col = net, group = net),
            linewidth = 1) +
  coord_cartesian(xlim = c(intervention_date_T - 27, intervention_date_T + 365*3),
                  ylim = c(0, 0.7)) +
  geom_vline(xintercept = intervention_date_T, linetype = 2, linewidth = 1, alpha = 0.5) +
  xlab("Year") +
  ylab("Malaria prevalence in children\naged 0.5 to 14 years old") +
  theme_classic() +
  ggtitle("Tanzania counterfactual") +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 14),
        legend.title = element_text(size=14),
        legend.background = element_blank(),
        legend.position = c(0.8, 0.15)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_x_date(expand = c(0, 0)) +
  scale_colour_manual(labels = c("No trial nets added",
                                 "Trial nets replaced\nwith pyrethroid-only nets",
                                 "Trial nets not replaced"
  ), 
  values = c("black", "darkgreen", "skyblue"),
  name = "",
  breaks = c("None", "IG2", "IG2_tu")) +
  scale_fill_manual(labels = c("No trial nets added",
                               "Trial nets replaced\nwith pyrethroid-only nets",
                               "Trial nets not replaced"
  ), 
  values = c("black", "darkgreen", "skyblue"),
  name = "",
  breaks = c("None", "IG2", "IG2_tu"))


cf_plot_B <- ggplot() +
  geom_ribbon(data = subset(B_tu_df, date <= intervention_date_B),
              aes(x = date, ymin = upper, ymax = lower, fill = net, group = net),
              alpha = 0.1) +
  
  geom_line(data = subset(B_tu_df, date <= intervention_date_B),
            aes(x = date, y = prev, col = net, group = net),
            alpha = 0.1, linewidth = 0.2) +
  
  geom_ribbon(data = subset(B_tu_df, date >= intervention_date_B),
              aes(x = date, ymin = upper, ymax = lower, fill = net, group = net),
              alpha = 0.25) +
  
  geom_line(data = subset(B_tu_df, date >= intervention_date_B),
            aes(x = date, y = prev, col = net, group = net),
            linewidth = 1) +
  
  coord_cartesian(xlim = c(intervention_date_B - 20, intervention_date_B + 365*3),
                  ylim = c(0, 0.7)) +
  geom_vline(xintercept = intervention_date_B, linetype = 2, linewidth = 1, alpha = 0.5) +
  xlab("Year") +
  ylab("Malaria prevalence in\npeople of all ages") +
  theme_classic() +
  ggtitle("Benin counterfactual") +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 14),
        legend.title = element_text(size=14),
        legend.background = element_blank(),
        legend.position = c(0.8, 0.85)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_x_date(expand = c(0, 0)) +
  scale_colour_manual(labels = c("No trial nets added",
                                 "Trial nets replaced\nwith pyrethroid-only nets",
                                 "Trial nets not replaced"
  ), 
  values = c("black", "darkgreen", "skyblue"),
  name = "",
  breaks = c("None", "IG2", "IG2_tu")) +
  scale_fill_manual(labels = c("No trial nets added",
                               "Trial nets replaced\nwith pyrethroid-only nets",
                               "Trial nets not replaced"
  ), 
  values = c("black", "darkgreen", "skyblue"),
  name = "",
  breaks = c("None", "IG2", "IG2_tu"))

return(list("T_tu_df" = T_tu_df,
            "B_tu_df" = B_tu_df,
            "cf_plot_T" = cf_plot_T,
            "cf_plot_B" = cf_plot_B))

}

#####################################################################
###### plotting the prevalence values for the trial simulations #####
#####################################################################

arm_sims_all <- extract_prev_arms(r_in = "o")

arm_sims_all_ll <- extract_prev_arms(r_in = "ll")

tu_sims_all <- extract_prev_diff_top_up(r_in = "o")

tu_sims_all_ll <- extract_prev_diff_top_up(r_in = "ll")

#################################
##### actual vs fitted plot #####
#################################

arm_sims_T_o <- arm_sims_all$arm_sims_T
arm_sims_B_o <- arm_sims_all$arm_sims_T
  
arm_sims_T_ll <- arm_sims_all_ll$arm_sims_T
arm_sims_B_ll <- arm_sims_all_ll$arm_sims_B

inds_T_o <- match(interaction(df_T$date, df_T$Net), interaction(arm_sims_T_o$date, arm_sims_T_o$net))
inds_T_ll <- match(interaction(df_T$date, df_T$Net), interaction(arm_sims_T_ll$date, arm_sims_T_ll$net))

inds_B_o <- match(interaction(df_B$date, df_B$Net), interaction(arm_sims_B_o$date, arm_sims_B_o$net))
inds_B_ll <- match(interaction(df_B$date, df_B$Net), interaction(arm_sims_B_ll$date, arm_sims_B_ll$net))
  
df_T_af <- as.data.frame(df_T) %>% mutate(pred_prev = arm_sims_T_ll[inds_T_ll, "prev"],
                                          pred_prev_l = arm_sims_T_o[inds_T_o, "lower"],
                                          pred_prev_u = arm_sims_T_o[inds_T_o, "upper"],
                                          pred_prev_l_ll = arm_sims_T_ll[inds_T_ll, "lower"],
                                          pred_prev_u_ll = arm_sims_T_ll[inds_T_ll, "upper"])

df_B_af <- as.data.frame(df_B) %>% mutate(pred_prev = arm_sims_B_ll[inds_B_ll, "prev"],
                                          pred_prev_l = arm_sims_B_o[inds_B_o, "lower"],
                                          pred_prev_u = arm_sims_B_o[inds_B_o, "upper"],
                                          pred_prev_l_ll = arm_sims_B_ll[inds_B_ll, "lower"],
                                          pred_prev_u_ll = arm_sims_B_ll[inds_B_ll, "upper"])
  
df_af <- rbind(df_T_af[, c("Net", "Country", "date", "Malaria_prevalence", "Malaria_prevalence_u", 
                        "Malaria_prevalence_l", "pred_prev", "pred_prev_l", "pred_prev_u", "pred_prev_l_ll", "pred_prev_u_ll")], 
               df_B_af[, c("Net", "Country", "date", "Malaria_prevalence", "Malaria_prevalence_u",
                        "Malaria_prevalence_l", "pred_prev", "pred_prev_l", "pred_prev_u", "pred_prev_l_ll", "pred_prev_u_ll")])

df_af <- df_af[is.na(df_af$Malaria_prevalence)!=1, ]
  
df_af_plot <- subset(df_af, !(date %in% c(baseline_start_date_B,
                                          baseline_start_date_T)) & Net != "RG") # excluding calibrations
  
af_plot <- ggplot(data = df_af_plot) +
  geom_abline(slope = 1, linetype = 2, linewidth = 1.1) +
  # geom_smooth(formula = y ~ x-1, se = TRUE, aes(x = Malaria_prevalence, y = pred_prev),
  #              method = "lm", col = "grey50", fullrange = TRUE, alpha = 0.1) +
  geom_pointrange(size = 1, 
                  aes(x = Malaria_prevalence, y = pred_prev, ymin = pred_prev_u, ymax = pred_prev_l, col = Net, shape = Country)) +
  geom_errorbarh(aes(y = pred_prev, xmin = Malaria_prevalence_l, xmax = Malaria_prevalence_u, col = Net, group = Country)) +
  geom_linerange(aes(x = Malaria_prevalence, ymin = pred_prev_u_ll, ymax = pred_prev_l_ll, col = Net, group = Country), linetype = 3, linewidth = 0.85) +
  theme_classic() + 
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 14),
        legend.title = element_text(size=14),
        legend.background = element_blank(),
        legend.position = c(0.85, 0.25),
        legend.box.background = element_rect(colour = "black")) +
  scale_colour_manual(values = c("blue", "aquamarine", "darkgreen"),
                      labels = c("pyrethroid-only", "pyrethroid-PBO", "pyrethroid-pyrrole"),
                      breaks =c("Pyrethroid_only", "PBO", "IG2")) +
  scale_shape_manual(values = c(15, 16)) +
  ylab("Predicted malaria prevalence") + xlab("Observed malaria prevalence") +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent, breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent, breaks = seq(0, 1, 0.1)) +
  coord_cartesian(xlim = c(0, 0.7), ylim = c(0, 0.7)) 


res_sq <- with(df_af_plot, (Malaria_prevalence - pred_prev)^2)
tot_sq <- with(df_af_plot, (Malaria_prevalence - mean(Malaria_prevalence))^2)
rs <- 1 - sum(res_sq)/sum(tot_sq)

af_plot + ggtitle(bquote(R^2==.(round(rs, digits = 3))))

#####################
##### incidence #####
#####################

calc_incidence <- function(output, s_date, e_date, min_age, max_age){
  s_time <- which(output$date == s_date)
  e_time <- which(output$date == e_date)
  return(output[s_time:e_time,paste0("n_inc_clinical_",min_age,"_",max_age)]/
           output[s_time:e_time, paste0("n_", min_age,"_", max_age)])
}

# returns the incidence per person year
extract_inc <- function(Location_in,
                        Net_in,
                        top_up_name,
                        Country_in,
                        r_model,
                        int_Net_in,
                        top_up_net,
                        trial_net_cov_in,
                        vals_pred, 
                        pred_prev,
                        n_rep = n_mcmc,
                        s_date, 
                        e_date,
                        min_age,
                        max_age){
  
  #p <- if(bioassay_uncertainty == "middle"){0.5} else if(bioassay_uncertainty == "lower"){0.025} else{0.975}
  
  inds <- which(vals_pred$Location_in == Location_in &
                  vals_pred$Net_in == Net_in &
                  vals_pred$top_up_name == top_up_name &
                  vals_pred$Country_in == Country_in &
                  vals_pred$r_model == r_model &
                  vals_pred$int_Net_in == int_Net_in &
                  vals_pred$top_up_net == top_up_net &
                  is.na(vals_pred$trial_net_cov_in) == is.na(trial_net_cov_in))
  
  if(length(inds) != n_rep){return(NA)}else{
    
    mat <- sapply(inds, function(i, p){calc_incidence(output = p[[i]],
                                                      s_date = s_date, 
                                                      e_date = e_date,
                                                      min_age = min_age, 
                                                      max_age = max_age)}, 
                  p = pred_prev)
  }
  
  
  
  if(nrow(mat)>1){
    mat <- colSums(mat)
  }
  
  mat <- mat / (difftime(e_date, s_date, units = c("days"))[[1]] / 365) # changing to the per person per year
  
  return(mat) 
}

# extracting the incidence values
calc_rel_reduction <- function(Location_in_trial,
                               Location_in_p,
                               Net_in,
                               top_up_name_trial,
                               top_up_name_p,
                               Country_in,
                               int_Net_in,
                               top_up_net,
                               trial_net_cov_in,
                               vals_pred, 
                               pred_prev,
                               n_rep = n_mcmc,
                               s_date = intervention_date_T, 
                               e_date = intervention_date_T + 365,
                               min_age = 0.5*365,
                               max_age = 10*365,
                               r_in,
                               quant){
  
  # p_only refers to pyrethroid only net in the pyrethroid-only arm
  t_inc <- extract_inc(Location_in = Location_in_trial,
                       Net_in = Net_in, 
                       top_up_name = top_up_name_trial,
                       Country_in = Country_in,
                       r_model = r_in,
                       int_Net_in = int_Net_in,
                       top_up_net = top_up_net,
                       trial_net_cov_in = trial_net_cov_in,
                       vals_pred = vals_pred, 
                       pred_prev = pred_prev,
                       n_rep = n_rep,
                       s_date = s_date, 
                       e_date = e_date,
                       min_age = min_age,
                       max_age = max_age)
  
  p_inc <- extract_inc(Location_in = Location_in_p,
                               Net_in = "Pyrethroid_only", 
                               top_up_name = top_up_name_p,
                               Country_in = Country_in,
                               r_model = r_in,
                               int_Net_in = "Pyrethroid_only",
                               top_up_net = "Pyrethroid_only",
                               trial_net_cov_in = NA,
                               vals_pred = vals_pred, 
                               pred_prev = pred_prev,
                               n_rep = n_rep,
                               s_date = s_date, 
                               e_date = e_date,
                               min_age = min_age,
                               max_age = max_age)
  
  out <- (p_inc - t_inc)/p_inc
  
  return(unname(quantile(out, probs = quant)))
}  

calc_rel_reduction_cf <- function(Location_in_trial,
                                  Net_in,
                                  top_up_name_trial,
                                  Country_in,
                                  int_Net_in,
                                  top_up_net,
                                  trial_net_cov_in,
                                  vals_pred, 
                                  pred_prev,
                                  n_rep = n_mcmc,
                                  s_date = intervention_date_T, 
                                  e_date = intervention_date_T + 365,
                                  min_age = 0.5*365,
                                  max_age = 10*365,
                                  r_in,
                                  quant){
  
  # p_only refers to pyrethroid only net in the pyrethroid-only arm
  t_inc <- extract_inc(Location_in = Location_in_trial,
                       Net_in = Net_in, 
                       top_up_name = top_up_name_trial,
                       Country_in = Country_in,
                       r_model = r_in,
                       int_Net_in = int_Net_in,
                       top_up_net = top_up_net,
                       trial_net_cov_in = trial_net_cov_in,
                       vals_pred = vals_pred, 
                       pred_prev = pred_prev,
                       n_rep = n_rep,
                       s_date = s_date, 
                       e_date = e_date,
                       min_age = min_age,
                       max_age = max_age)
  
  p_inc <- extract_inc(Location_in = Location_in_trial,
                       Net_in = Net_in, 
                       top_up_name = top_up_name_trial,
                       Country_in = Country_in,
                       r_model = r_in,
                       int_Net_in = int_Net_in,
                       top_up_net = "none",
                       trial_net_cov_in = 0,
                       vals_pred = vals_pred, 
                       pred_prev = pred_prev,
                       n_rep = n_rep,
                       s_date = s_date, 
                       e_date = e_date,
                       min_age = min_age,
                       max_age = max_age)
  
  out <- (p_inc - t_inc)/p_inc
  
  return(unname(quantile(out, probs = quant)))
}  

rel_reduction_T <- data.frame(s_times = rep(c(intervention_date_T, 
                                                intervention_date_T + 365,
                                                intervention_date_T + 365 * 2,
                                                intervention_date_T,
                                                intervention_date_T),
                                                3),
                                  e_times = rep(c(intervention_date_T + 365,
                                                  intervention_date_T + 365 * 2,
                                                  intervention_date_T + 365 * 3,
                                                  intervention_date_T + 365 * 2,
                                                  intervention_date_T + 365 * 3),
                                                3),
                                  year = rep(c("Year 1", "Year 2", "Year 3", "Years 1 & 2", "Years 1, 2 & 3"), 3),
                                  net = sort(rep(c("Pyrethroid_only", "PBO", "IG2"), 5))
                              ) %>% rowwise() %>%
  
  mutate(inc_ll_m = extract_inc(Location_in = "Misungwi",
                           Net_in = net, 
                           top_up_name = ifelse(net == "IG2", "top_up_IG2_T", ifelse(net == "PBO", "top_up_PBO_T", "top_up_p_only_T")),
                           Country_in = "Tanzania",
                           r_model = "ll",
                           int_Net_in = net,
                           top_up_net = "Pyrethroid_only",
                           trial_net_cov_in = NA,
                           vals_pred = T_vals_pred, 
                           pred_prev = pred_prev_T,
                           n_rep = n_mcmc,
                           s_date = s_times, 
                           e_date = e_times,
                           min_age = 0.5*365,
                           max_age = 10*365) %>% quantile(probs = 0.5) %>% unname(),
         
         inc_ll_l = extract_inc(Location_in = "Misungwi",
                                Net_in = net, 
                                top_up_name = ifelse(net == "IG2", "top_up_IG2_T", ifelse(net == "PBO", "top_up_PBO_T", "top_up_p_only_T")),
                                Country_in = "Tanzania",
                                r_model = "ll",
                                int_Net_in = net,
                                top_up_net = "Pyrethroid_only",
                                trial_net_cov_in = NA,
                                vals_pred = T_vals_pred, 
                                pred_prev = pred_prev_T,
                                n_rep = n_mcmc,
                                s_date = s_times, 
                                e_date = e_times,
                                min_age = 0.5*365,
                                max_age = 10*365) %>% quantile(probs = 0.025) %>% unname(),
         
         inc_ll_u = extract_inc(Location_in = "Misungwi",
                                Net_in = net, 
                                top_up_name = ifelse(net == "IG2", "top_up_IG2_T", ifelse(net == "PBO", "top_up_PBO_T", "top_up_p_only_T")),
                                Country_in = "Tanzania",
                                r_model = "ll",
                                int_Net_in = net,
                                top_up_net = "Pyrethroid_only",
                                trial_net_cov_in = NA,
                                vals_pred = T_vals_pred, 
                                pred_prev = pred_prev_T,
                                n_rep = n_mcmc,
                                s_date = s_times, 
                                e_date = e_times,
                                min_age = 0.5*365,
                                max_age = 10*365) %>% quantile(probs = 0.975) %>% unname(),
         
         tu_inc_ll_m = extract_inc(Location_in = "Misungwi",
                                Net_in = net, 
                                top_up_name = ifelse(net == "IG2", "top_up_IG2_T", ifelse(net == "PBO", "top_up_PBO_T", "top_up_p_only_T")),
                                Country_in = "Tanzania",
                                r_model = "ll",
                                int_Net_in = net,
                                top_up_net = net,
                                trial_net_cov_in = NA,
                                vals_pred = T_vals_pred, 
                                pred_prev = pred_prev_T,
                                n_rep = n_mcmc,
                                s_date = s_times, 
                                e_date = e_times,
                                min_age = 0.5*365,
                                max_age = 10*365) %>% quantile(probs = 0.5) %>% unname(),
         
         tu_inc_ll_l = extract_inc(Location_in = "Misungwi",
                                Net_in = net, 
                                top_up_name = ifelse(net == "IG2", "top_up_IG2_T", ifelse(net == "PBO", "top_up_PBO_T", "top_up_p_only_T")),
                                Country_in = "Tanzania",
                                r_model = "ll",
                                int_Net_in = net,
                                top_up_net = net,
                                trial_net_cov_in = NA,
                                vals_pred = T_vals_pred, 
                                pred_prev = pred_prev_T,
                                n_rep = n_mcmc,
                                s_date = s_times, 
                                e_date = e_times,
                                min_age = 0.5*365,
                                max_age = 10*365) %>% quantile(probs = 0.025) %>% unname(),
         
         tu_inc_ll_u = extract_inc(Location_in = "Misungwi",
                                Net_in = net, 
                                top_up_name = ifelse(net == "IG2", "top_up_IG2_T", ifelse(net == "PBO", "top_up_PBO_T", "top_up_p_only_T")),
                                Country_in = "Tanzania",
                                r_model = "ll",
                                int_Net_in = net,
                                top_up_net = net,
                                trial_net_cov_in = NA,
                                vals_pred = T_vals_pred, 
                                pred_prev = pred_prev_T,
                                n_rep = n_mcmc,
                                s_date = s_times, 
                                e_date = e_times,
                                min_age = 0.5*365,
                                max_age = 10*365) %>% quantile(probs = 0.975) %>% unname(),
         
         rr_m = calc_rel_reduction(Location_in_trial = "Misungwi",
                            Location_in_p = "Misungwi",
                            Net_in = net,
                            top_up_name_trial = ifelse(net == "IG2", "top_up_IG2_T", ifelse(net == "PBO", "top_up_PBO_T", "top_up_p_only_T")),
                            top_up_name_p = "top_up_p_only_T",
                            Country_in = "Tanzania",
                            int_Net_in = net,
                            top_up_net = "Pyrethroid_only",
                            trial_net_cov_in = NA,
                            vals_pred = T_vals_pred, 
                            pred_prev = pred_prev_T,
                            n_rep = n_mcmc,
                            s_date = s_times, 
                            e_date = e_times,
                            min_age = 0.5*365,
                            max_age = 10*365,
                            r_in = "ll",
                            quant = 0.5),
         
         rr_l = calc_rel_reduction(Location_in_trial = "Misungwi",
                                   Location_in_p = "Misungwi",
                                   Net_in = net,
                                   top_up_name_trial = ifelse(net == "IG2", "top_up_IG2_T", ifelse(net == "PBO", "top_up_PBO_T", "top_up_p_only_T")),
                                   top_up_name_p = "top_up_p_only_T",
                                   Country_in = "Tanzania",
                                   int_Net_in = net,
                                   top_up_net = "Pyrethroid_only",
                                   trial_net_cov_in = NA,
                                   vals_pred = T_vals_pred, 
                                   pred_prev = pred_prev_T,
                                   n_rep = n_mcmc,
                                   s_date = s_times, 
                                   e_date = e_times,
                                   min_age = 0.5*365,
                                   max_age = 10*365,
                                   r_in = "ll",
                                   quant = 0.025),
         
         rr_u = calc_rel_reduction(Location_in_trial = "Misungwi",
                                   Location_in_p = "Misungwi",
                                   Net_in = net,
                                   top_up_name_trial = ifelse(net == "IG2", "top_up_IG2_T", ifelse(net == "PBO", "top_up_PBO_T", "top_up_p_only_T")),
                                   top_up_name_p = "top_up_p_only_T",
                                   Country_in = "Tanzania",
                                   int_Net_in = net,
                                   top_up_net = "Pyrethroid_only",
                                   trial_net_cov_in = NA,
                                   vals_pred = T_vals_pred, 
                                   pred_prev = pred_prev_T,
                                   n_rep = n_mcmc,
                                   s_date = s_times, 
                                   e_date = e_times,
                                   min_age = 0.5*365,
                                   max_age = 10*365,
                                   r_in = "ll",
                                   quant = 0.975),
         
         tu_rr_m = calc_rel_reduction(Location_in_trial = "Misungwi",
                                   Location_in_p = "Misungwi",
                                   Net_in = net,
                                   top_up_name_trial = ifelse(net == "IG2", "top_up_IG2_T", ifelse(net == "PBO", "top_up_PBO_T", "top_up_p_only_T")),
                                   top_up_name_p = "top_up_p_only_T",
                                   Country_in = "Tanzania",
                                   int_Net_in = net,
                                   top_up_net = net,
                                   trial_net_cov_in = NA,
                                   vals_pred = T_vals_pred, 
                                   pred_prev = pred_prev_T,
                                   n_rep = n_mcmc,
                                   s_date = s_times, 
                                   e_date = e_times,
                                   min_age = 0.5*365,
                                   max_age = 10*365,
                                   r_in = "ll",
                                   quant = 0.5),
         
         tu_rr_l = calc_rel_reduction(Location_in_trial = "Misungwi",
                                   Location_in_p = "Misungwi",
                                   Net_in = net,
                                   top_up_name_trial = ifelse(net == "IG2", "top_up_IG2_T", ifelse(net == "PBO", "top_up_PBO_T", "top_up_p_only_T")),
                                   top_up_name_p = "top_up_p_only_T",
                                   Country_in = "Tanzania",
                                   int_Net_in = net,
                                   top_up_net = net,
                                   trial_net_cov_in = NA,
                                   vals_pred = T_vals_pred, 
                                   pred_prev = pred_prev_T,
                                   n_rep = n_mcmc,
                                   s_date = s_times, 
                                   e_date = e_times,
                                   min_age = 0.5*365,
                                   max_age = 10*365,
                                   r_in = "ll",
                                   quant = 0.025),
         
         tu_rr_u = calc_rel_reduction(Location_in_trial = "Misungwi",
                                   Location_in_p = "Misungwi",
                                   Net_in = net,
                                   top_up_name_trial = ifelse(net == "IG2", "top_up_IG2_T", ifelse(net == "PBO", "top_up_PBO_T", "top_up_p_only_T")),
                                   top_up_name_p = "top_up_p_only_T",
                                   Country_in = "Tanzania",
                                   int_Net_in = net,
                                   top_up_net = net,
                                   trial_net_cov_in = NA,
                                   vals_pred = T_vals_pred, 
                                   pred_prev = pred_prev_T,
                                   n_rep = n_mcmc,
                                   s_date = s_times, 
                                   e_date = e_times,
                                   min_age = 0.5*365,
                                   max_age = 10*365,
                                   r_in = "ll",
                                   quant = 0.975),
         
         overall_efficacy_m = calc_rel_reduction_cf(Location_in_trial = "Misungwi",
                                                    Net_in = net,
                                                    top_up_name_trial = ifelse(net == "IG2", "top_up_IG2_T", ifelse(net == "PBO", "top_up_PBO_T", "top_up_p_only_T")),
                                                    Country_in = "Tanzania",
                                                    int_Net_in = net,
                                                    top_up_net = "Pyrethroid_only",
                                                    trial_net_cov_in = NA,
                                                    vals_pred = T_vals_pred, 
                                                    pred_prev = pred_prev_T,
                                                    n_rep = n_mcmc,
                                                    s_date = s_times, 
                                                    e_date = e_times,
                                                    min_age = 0.5*365,
                                                    max_age = 10*365,
                                                    r_in = "ll",
                                                    quant = 0.5),
         
         overall_efficacy_l = calc_rel_reduction_cf(Location_in_trial = "Misungwi",
                                                    Net_in = net,
                                                    top_up_name_trial = ifelse(net == "IG2", "top_up_IG2_T", ifelse(net == "PBO", "top_up_PBO_T", "top_up_p_only_T")),
                                                    Country_in = "Tanzania",
                                                    int_Net_in = net,
                                                    top_up_net = "Pyrethroid_only",
                                                    trial_net_cov_in = NA,
                                                    vals_pred = T_vals_pred, 
                                                    pred_prev = pred_prev_T,
                                                    n_rep = n_mcmc,
                                                    s_date = s_times, 
                                                    e_date = e_times,
                                                    min_age = 0.5*365,
                                                    max_age = 10*365,
                                                    r_in = "ll",
                                                    quant = 0.025),
         
         overall_efficacy_u = calc_rel_reduction_cf(Location_in_trial = "Misungwi",
                                                    Net_in = net,
                                                    top_up_name_trial = ifelse(net == "IG2", "top_up_IG2_T", ifelse(net == "PBO", "top_up_PBO_T", "top_up_p_only_T")),
                                                    Country_in = "Tanzania",
                                                    int_Net_in = net,
                                                    top_up_net = "Pyrethroid_only",
                                                    trial_net_cov_in = NA,
                                                    vals_pred = T_vals_pred, 
                                                    pred_prev = pred_prev_T,
                                                    n_rep = n_mcmc,
                                                    s_date = s_times, 
                                                    e_date = e_times,
                                                    min_age = 0.5*365,
                                                    max_age = 10*365,
                                                    r_in = "ll",
                                                    quant = 0.975),
         
         tu_overall_efficacy_m = calc_rel_reduction_cf(Location_in_trial = "Misungwi",
                                                      Net_in = net,
                                                      top_up_name_trial = ifelse(net == "IG2", "top_up_IG2_T", ifelse(net == "PBO", "top_up_PBO_T", "top_up_p_only_T")),
                                                      Country_in = "Tanzania",
                                                      int_Net_in = net,
                                                      top_up_net = net,
                                                      trial_net_cov_in = NA,
                                                      vals_pred = T_vals_pred, 
                                                      pred_prev = pred_prev_T,
                                                      n_rep = n_mcmc,
                                                      s_date = s_times, 
                                                      e_date = e_times,
                                                      min_age = 0.5*365,
                                                      max_age = 10*365,
                                                      r_in = "ll",
                                                      quant = 0.5),
         
         tu_overall_efficacy_l = calc_rel_reduction_cf(Location_in_trial = "Misungwi",
                                                    Net_in = net,
                                                    top_up_name_trial = ifelse(net == "IG2", "top_up_IG2_T", ifelse(net == "PBO", "top_up_PBO_T", "top_up_p_only_T")),
                                                    Country_in = "Tanzania",
                                                    int_Net_in = net,
                                                    top_up_net = net,
                                                    trial_net_cov_in = NA,
                                                    vals_pred = T_vals_pred, 
                                                    pred_prev = pred_prev_T,
                                                    n_rep = n_mcmc,
                                                    s_date = s_times, 
                                                    e_date = e_times,
                                                    min_age = 0.5*365,
                                                    max_age = 10*365,
                                                    r_in = "ll",
                                                    quant = 0.025),
         
         tu_overall_efficacy_u = calc_rel_reduction_cf(Location_in_trial = "Misungwi",
                                                    Net_in = net,
                                                    top_up_name_trial = ifelse(net == "IG2", "top_up_IG2_T", ifelse(net == "PBO", "top_up_PBO_T", "top_up_p_only_T")),
                                                    Country_in = "Tanzania",
                                                    int_Net_in = net,
                                                    top_up_net = net,
                                                    trial_net_cov_in = NA,
                                                    vals_pred = T_vals_pred, 
                                                    pred_prev = pred_prev_T,
                                                    n_rep = n_mcmc,
                                                    s_date = s_times, 
                                                    e_date = e_times,
                                                    min_age = 0.5*365,
                                                    max_age = 10*365,
                                                    r_in = "ll",
                                                    quant = 0.975),
         Country = "Tanzania")



  
rel_reduction_B <- data.frame(s_times = rep(c(intervention_date_B, 
                           intervention_date_B + 365,
                           intervention_date_B + 365 * 2,
                           intervention_date_B,
                           intervention_date_B), 2),
           e_times = rep(c(intervention_date_B + 365,
                           intervention_date_B + 365 * 2,
                           intervention_date_B + 365 * 3,
                           intervention_date_B + 365 * 2,
                           intervention_date_B + 365 * 3), 2),
           year = rep(c("Year 1", "Year 2", "Year 3", "Years 1 & 2", "Years 1, 2 & 3"), 2),
           net = sort(rep(c("Pyrethroid_only", "IG2"), 5))) %>% 
  rowwise() %>% 
  mutate(inc_ll_m = extract_inc(Location_in = ifelse(net == "IG2", "Zagnanado", "Ouinhi"),
                                Net_in = net, 
                                top_up_name = ifelse(net == "IG2", "top_up_IG2_B", "top_up_p_only_B"),
                                Country_in = "Benin",
                                r_model = "ll",
                                int_Net_in = net,
                                top_up_net = "Pyrethroid_only",
                                trial_net_cov_in = NA,
                                vals_pred = B_vals_pred, 
                                pred_prev = pred_prev_B,
                                n_rep = n_mcmc,
                                s_date = s_times, 
                                e_date = e_times,
                                min_age = 0.5*365,
                                max_age = 10*365) %>% quantile(probs = 0.5) %>% unname(),
         
         inc_ll_l = extract_inc(Location_in = ifelse(net == "IG2", "Zagnanado", "Ouinhi"),
                                Net_in = net, 
                                top_up_name = ifelse(net == "IG2", "top_up_IG2_B", "top_up_p_only_B"),
                                Country_in = "Benin",
                                r_model = "ll",
                                int_Net_in = net,
                                top_up_net = "Pyrethroid_only",
                                trial_net_cov_in = NA,
                                vals_pred = B_vals_pred, 
                                pred_prev = pred_prev_B,
                                n_rep = n_mcmc,
                                s_date = s_times, 
                                e_date = e_times,
                                min_age = 0.5*365,
                                max_age = 10*365) %>% quantile(probs = 0.025) %>% unname(),
         
         inc_ll_u = extract_inc(Location_in = ifelse(net == "IG2", "Zagnanado", "Ouinhi"),
                                Net_in = net, 
                                top_up_name = ifelse(net == "IG2", "top_up_IG2_B", "top_up_p_only_B"),
                                Country_in = "Benin",
                                r_model = "ll",
                                int_Net_in = net,
                                top_up_net = "Pyrethroid_only",
                                trial_net_cov_in = NA,
                                vals_pred = B_vals_pred, 
                                pred_prev = pred_prev_B,
                                n_rep = n_mcmc,
                                s_date = s_times, 
                                e_date = e_times,
                                min_age = 0.5*365,
                                max_age = 10*365) %>% quantile(probs = 0.975) %>% unname(),
         
         tu_inc_ll_m = extract_inc(Location_in = ifelse(net == "IG2", "Zagnanado", "Ouinhi"),
                                   Net_in = net, 
                                   top_up_name = ifelse(net == "IG2", "top_up_IG2_B", "top_up_p_only_B"),
                                   Country_in = "Benin",
                                   r_model = "ll",
                                   int_Net_in = net,
                                   top_up_net = net,
                                   trial_net_cov_in = NA,
                                   vals_pred = B_vals_pred, 
                                   pred_prev = pred_prev_B,
                                   n_rep = n_mcmc,
                                   s_date = s_times, 
                                   e_date = e_times,
                                   min_age = 0.5*365,
                                   max_age = 10*365) %>% quantile(probs = 0.5) %>% unname(),
         
         tu_inc_ll_l = extract_inc(Location_in = ifelse(net == "IG2", "Zagnanado", "Ouinhi"),
                                   Net_in = net, 
                                   top_up_name = ifelse(net == "IG2", "top_up_IG2_B", "top_up_p_only_B"),
                                   Country_in = "Benin",
                                   r_model = "ll",
                                   int_Net_in = net,
                                   top_up_net = net,
                                   trial_net_cov_in = NA,
                                   vals_pred = B_vals_pred, 
                                   pred_prev = pred_prev_B,
                                   n_rep = n_mcmc,
                                   s_date = s_times, 
                                   e_date = e_times,
                                   min_age = 0.5*365,
                                   max_age = 10*365) %>% quantile(probs = 0.025) %>% unname(),
         
         tu_inc_ll_u = extract_inc(Location_in = ifelse(net == "IG2", "Zagnanado", "Ouinhi"),
                                   Net_in = net, 
                                   top_up_name = ifelse(net == "IG2", "top_up_IG2_B", "top_up_p_only_B"),
                                   Country_in = "Benin",
                                   r_model = "ll",
                                   int_Net_in = net,
                                   top_up_net = net,
                                   trial_net_cov_in = NA,
                                   vals_pred = B_vals_pred, 
                                   pred_prev = pred_prev_B,
                                   n_rep = n_mcmc,
                                   s_date = s_times, 
                                   e_date = e_times,
                                   min_age = 0.5*365,
                                   max_age = 10*365) %>% quantile(probs = 0.975) %>% unname(),
         
         rr_m = calc_rel_reduction(Location_in_trial = ifelse(net == "IG2", "Zagnanado", "Ouinhi"),
                                   Location_in_p = "Ouinhi",
                                   Net_in = net,
                                   top_up_name_trial = ifelse(net == "IG2", "top_up_IG2_B", "top_up_p_only_B"),
                                   top_up_name_p = "top_up_p_only_B",
                                   Country_in = "Benin",
                                   int_Net_in = net,
                                   top_up_net = "Pyrethroid_only",
                                   trial_net_cov_in = NA,
                                   vals_pred = B_vals_pred, 
                                   pred_prev = pred_prev_B,
                                   n_rep = n_mcmc,
                                   s_date = s_times, 
                                   e_date = e_times,
                                   min_age = 0.5*365,
                                   max_age = 10*365,
                                   r_in = "ll",
                                   quant = 0.5),
         
         rr_l = calc_rel_reduction(Location_in_trial = ifelse(net == "IG2", "Zagnanado", "Ouinhi"),
                                   Location_in_p = "Ouinhi",
                                   Net_in = net,
                                   top_up_name_trial = ifelse(net == "IG2", "top_up_IG2_B", "top_up_p_only_B"),
                                   top_up_name_p = "top_up_p_only_B",
                                   Country_in = "Benin",
                                   int_Net_in = net,
                                   top_up_net = "Pyrethroid_only",
                                   trial_net_cov_in = NA,
                                   vals_pred = B_vals_pred, 
                                   pred_prev = pred_prev_B,
                                   n_rep = n_mcmc,
                                   s_date = s_times, 
                                   e_date = e_times,
                                   min_age = 0.5*365,
                                   max_age = 10*365,
                                   r_in = "ll",
                                   quant = 0.025),
         
         rr_u = calc_rel_reduction(Location_in_trial = ifelse(net == "IG2", "Zagnanado", "Ouinhi"),
                                   Location_in_p = "Ouinhi",
                                   Net_in = net,
                                   top_up_name_trial = ifelse(net == "IG2", "top_up_IG2_B", "top_up_p_only_B"),
                                   top_up_name_p = "top_up_p_only_B",
                                   Country_in = "Benin",
                                   int_Net_in = net,
                                   top_up_net = "Pyrethroid_only",
                                   trial_net_cov_in = NA,
                                   vals_pred = B_vals_pred, 
                                   pred_prev = pred_prev_B,
                                   n_rep = n_mcmc,
                                   s_date = s_times, 
                                   e_date = e_times,
                                   min_age = 0.5*365,
                                   max_age = 10*365,
                                   r_in = "ll",
                                   quant = 0.975),
         
         tu_rr_m = calc_rel_reduction(Location_in_trial = ifelse(net == "IG2", "Zagnanado", "Ouinhi"),
                                      Location_in_p = "Ouinhi",
                                      Net_in = net,
                                      top_up_name_trial = ifelse(net == "IG2", "top_up_IG2_B", "top_up_p_only_B"),
                                      top_up_name_p = "top_up_p_only_B",
                                      Country_in = "Benin",
                                      int_Net_in = net,
                                      top_up_net = net,
                                      trial_net_cov_in = NA,
                                      vals_pred = B_vals_pred, 
                                      pred_prev = pred_prev_B,
                                      n_rep = n_mcmc,
                                      s_date = s_times, 
                                      e_date = e_times,
                                      min_age = 0.5*365,
                                      max_age = 10*365,
                                      r_in = "ll",
                                      quant = 0.5),
         
         tu_rr_l = calc_rel_reduction(Location_in_trial = ifelse(net == "IG2", "Zagnanado", "Ouinhi"),
                                      Location_in_p = "Ouinhi",
                                      Net_in = net,
                                      top_up_name_trial = ifelse(net == "IG2", "top_up_IG2_B", "top_up_p_only_B"),
                                      top_up_name_p = "top_up_p_only_B",
                                      Country_in = "Benin",
                                      int_Net_in = net,
                                      top_up_net = net,
                                      trial_net_cov_in = NA,
                                      vals_pred = B_vals_pred, 
                                      pred_prev = pred_prev_B,
                                      n_rep = n_mcmc,
                                      s_date = s_times, 
                                      e_date = e_times,
                                      min_age = 0.5*365,
                                      max_age = 10*365,
                                      r_in = "ll",
                                      quant = 0.025),
         
         tu_rr_u = calc_rel_reduction(Location_in_trial = ifelse(net == "IG2", "Zagnanado", "Ouinhi"),
                                      Location_in_p = "Ouinhi",
                                      Net_in = net,
                                      top_up_name_trial = ifelse(net == "IG2", "top_up_IG2_B", "top_up_p_only_B"),
                                      top_up_name_p = "top_up_p_only_B",
                                      Country_in = "Benin",
                                      int_Net_in = net,
                                      top_up_net = net,
                                      trial_net_cov_in = NA,
                                      vals_pred = B_vals_pred, 
                                      pred_prev = pred_prev_B,
                                      n_rep = n_mcmc,
                                      s_date = s_times, 
                                      e_date = e_times,
                                      min_age = 0.5*365,
                                      max_age = 10*365,
                                      r_in = "ll",
                                      quant = 0.975),
         
         overall_efficacy_m = calc_rel_reduction_cf(Location_in_trial = ifelse(net == "IG2", "Zagnanado", "Ouinhi"),
                                                    Net_in = net,
                                                    top_up_name_trial = ifelse(net == "IG2", "top_up_IG2_B", "top_up_p_only_B"),
                                                    Country_in = "Benin",
                                                    int_Net_in = net,
                                                    top_up_net = "Pyrethroid_only",
                                                    trial_net_cov_in = NA,
                                                    vals_pred = B_vals_pred, 
                                                    pred_prev = pred_prev_B,
                                                    n_rep = n_mcmc,
                                                    s_date = s_times, 
                                                    e_date = e_times,
                                                    min_age = 0.5*365,
                                                    max_age = 10*365,
                                                    r_in = "ll",
                                                    quant = 0.5),
         
         overall_efficacy_l = calc_rel_reduction_cf(Location_in_trial = ifelse(net == "IG2", "Zagnanado", "Ouinhi"),
                                                    Net_in = net,
                                                    top_up_name_trial = ifelse(net == "IG2", "top_up_IG2_B", "top_up_p_only_B"),
                                                    Country_in = "Benin",
                                                    int_Net_in = net,
                                                    top_up_net = "Pyrethroid_only",
                                                    trial_net_cov_in = NA,
                                                    vals_pred = B_vals_pred, 
                                                    pred_prev = pred_prev_B,
                                                    n_rep = n_mcmc,
                                                    s_date = s_times, 
                                                    e_date = e_times,
                                                    min_age = 0.5*365,
                                                    max_age = 10*365,
                                                    r_in = "ll",
                                                    quant = 0.025),
         
         overall_efficacy_u = calc_rel_reduction_cf(Location_in_trial = ifelse(net == "IG2", "Zagnanado", "Ouinhi"),
                                                    Net_in = net,
                                                    top_up_name_trial = ifelse(net == "IG2", "top_up_IG2_B", "top_up_p_only_B"),
                                                    Country_in = "Benin",
                                                    int_Net_in = net,
                                                    top_up_net = "Pyrethroid_only",
                                                    trial_net_cov_in = NA,
                                                    vals_pred = B_vals_pred, 
                                                    pred_prev = pred_prev_B,
                                                    n_rep = n_mcmc,
                                                    s_date = s_times, 
                                                    e_date = e_times,
                                                    min_age = 0.5*365,
                                                    max_age = 10*365,
                                                    r_in = "ll",
                                                    quant = 0.975),
         
         tu_overall_efficacy_m = calc_rel_reduction_cf(Location_in_trial = ifelse(net == "IG2", "Zagnanado", "Ouinhi"),
                                                       Net_in = net,
                                                       top_up_name_trial = ifelse(net == "IG2", "top_up_IG2_B", "top_up_p_only_B"),
                                                       Country_in = "Benin",
                                                       int_Net_in = net,
                                                       top_up_net = net,
                                                       trial_net_cov_in = NA,
                                                       vals_pred = B_vals_pred, 
                                                       pred_prev = pred_prev_B,
                                                       n_rep = n_mcmc,
                                                       s_date = s_times, 
                                                       e_date = e_times,
                                                       min_age = 0.5*365,
                                                       max_age = 10*365,
                                                       r_in = "ll",
                                                       quant = 0.5),
         
         tu_overall_efficacy_l = calc_rel_reduction_cf(Location_in_trial = ifelse(net == "IG2", "Zagnanado", "Ouinhi"),
                                                       Net_in = net,
                                                       top_up_name_trial = ifelse(net == "IG2", "top_up_IG2_B", "top_up_p_only_B"),
                                                       Country_in = "Benin",
                                                       int_Net_in = net,
                                                       top_up_net = net,
                                                       trial_net_cov_in = NA,
                                                       vals_pred = B_vals_pred, 
                                                       pred_prev = pred_prev_B,
                                                       n_rep = n_mcmc,
                                                       s_date = s_times, 
                                                       e_date = e_times,
                                                       min_age = 0.5*365,
                                                       max_age = 10*365,
                                                       r_in = "ll",
                                                       quant = 0.025),
         
         tu_overall_efficacy_u = calc_rel_reduction_cf(Location_in_trial = ifelse(net == "IG2", "Zagnanado", "Ouinhi"),
                                                       Net_in = net,
                                                       top_up_name_trial = ifelse(net == "IG2", "top_up_IG2_B", "top_up_p_only_B"),
                                                       Country_in = "Benin",
                                                       int_Net_in = net,
                                                       top_up_net = net,
                                                       trial_net_cov_in = NA,
                                                       vals_pred = B_vals_pred, 
                                                       pred_prev = pred_prev_B,
                                                       n_rep = n_mcmc,
                                                       s_date = s_times, 
                                                       e_date = e_times,
                                                       min_age = 0.5*365,
                                                       max_age = 10*365,
                                                       r_in = "ll",
                                                       quant = 0.975),
         Country = "Benin") 

rel_reduction <- rbind(rel_reduction_B, rel_reduction_T)
  
rel_reduction$Country <- factor(rel_reduction$Country, levels=c("Tanzania", "Benin"))

actual_inc <- read_xlsx("data/actual_incidence_estimates.xlsx") %>% 
  mutate(rr_m = (Pyrethroid_only_incidence - Trial_incidence)/Pyrethroid_only_incidence,
                                                                           sample = "observed",
         rr_l = NA, rr_u = NA) 

actual_inc$Country <- factor(actual_inc$Country, levels=c("Tanzania", "Benin"))
  
# rel_reduction <- rbind(rel_reduction_T, rel_reduction_B) %>% rowwise() %>% 
#     mutate(rr_m = (p_only_m - t_net_m)/p_only_m,
#            rr_l = (p_only_l - t_net_l)/p_only_l,
#            rr_u = (p_only_u - t_net_u)/p_only_u,
#                            
#                            overall_efficacy_m = (cf_net_m - t_net_m) / cf_net_m,
#                            overall_efficacy_l = (cf_net_l - t_net_l) / cf_net_l,
#                            overall_efficacy_u = (cf_net_u - t_net_u) / cf_net_u,
#                            
#                            rr_tu_m = (p_only_m - tu_net_m)/p_only_m,
#                            rr_tu_l = (p_only_l - tu_net_l)/p_only_l,
#                            rr_tu_u = (p_only_u - tu_net_u)/p_only_u,
#                            
#                            overall_efficacy_tu_m = (cf_net_m - tu_net_m) / cf_net_m,
#                            overall_efficacy_tu_l = (cf_net_l - tu_net_l) / cf_net_l,
#                            overall_efficacy_tu_u = (cf_net_u - tu_net_u) / cf_net_u)

rr_plot_df <- rbind(subset(rel_reduction, Country == "Tanzania" & net != "Pyrethroid_only" & net!= "RG" & year != "Year 3" & year != "Years 1, 2 & 3" |
                                           Country == "Benin" & net == "IG2")[,c("Country", "year", "rr_m", "rr_l", "rr_u", "net")] %>% 
                                 mutate(sample = "predicted"),
                               subset(actual_inc[,c("Country", "year", "rr_m", "rr_l", "rr_u", "net", "sample")], Country == "Tanzania" & net != "Pyrethroid_only" & net!= "RG" | Country == "Benin" & net == "IG2")) %>% 
                    rowwise() %>%  
                    mutate(net = if(net == "IG2"){"pyrethroid-pyrrole"} else{"pyrethroid-PBO"},
                           year = str_replace(year, "Year ", ""),
                           year = str_replace(year, "Years ", ""))

rr_plot_df$year <- factor(rr_plot_df$year, levels = c("1", "2", "3", "1 & 2", "1, 2 & 3"))

rr_plot_ll <- ggplot(data = rr_plot_df, 
                  aes(x = year, y = rr_m, ymin = rr_l, ymax = rr_u, fill = net, group = interaction(net, sample))) + 
  ggpattern::geom_col_pattern(alpha = 0.4, 
                               position = position_dodge(width = 0.9),
                               pattern_colour = "black",
                               col = "black",
                               aes(pattern = sample)) +
  #geom_col(alpha = 0.4, position = position_dodge()) +
  geom_errorbar(position=position_dodge(width = 0.9)) +
  scale_pattern_manual(values = c(predicted = "circle", observed = "none"), name = "Sample") +
  facet_wrap(~Country + net, scales = "free_x") +
  xlab("Year post trial net distribution") + ylab("Efficacy (clinical incidence per person-year averted\nrelative to the value in the pyrethroid-only arm)") +
  theme_classic() + theme(text = element_text(size = 20),
                          axis.text.y = element_text(size = 17),
                          axis.text.x = element_text(size = 12.75),
                          legend.text = element_text(size = 14),
                          legend.title = element_text(size=14),
                          legend.background = element_blank()) +
  scale_y_continuous(labels = scales::percent) +
  scale_colour_manual(values = c("aquamarine", "darkgreen"), name = "Net") +
  scale_fill_manual(values = c("aquamarine", "darkgreen"), name = "Net") +
  theme(legend.position = c(0.175, 0.825))


png(file = "figures/model_simulations_beta_binomial.png", height = 1100, width = 2175)

(arm_sims_all_ll$T_arm_plot | arm_sims_all_ll$B_arm_plot | af_plot) / (rr_plot_ll |  tu_sims_all_ll$cf_plot_T | tu_sims_all_ll$cf_plot_B) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold'))

dev.off()

arm_sims_all$T_arm_plot | arm_sims_all$B_arm_plot

############################################
##### extracting the prevalence values #####
############################################

# extracting the values
inc_out <- 
  rel_reduction %>% rowwise() %>%
  mutate(observed_arm = paste0(round(inc_ll_m, digits = 2), " (", round(rr_m, digits = 2)*100,"%)"),
         overall_arm = paste0(round(overall_efficacy_m, digits = 2)*100, "%"),
         observed_tu = paste0(round(tu_inc_ll_m, digits = 2), " (", round(tu_rr_m, digits = 2)*100,"%)"),
         overall_tu = paste0(round(tu_overall_efficacy_m, digits = 2)*100, "%")) %>% 
  select(c(Country, year, net, observed_arm, overall_arm, observed_tu, overall_tu))

write.csv(inc_out, file = "data/inc_estimates.csv")

write.csv(actual_inc %>% mutate(inc_out = paste0(Trial_incidence, " (", round(rr_m, digits = 2)*100, "%)")), file = "data/format_actual_inc_estimates.csv")

extract_efficacy <- function(Location_in,
                             Location_in_p,
                             Net_in,
                         top_up_name,
                         top_up_name_p,
                         Country_in,
                         r_model,
                         int_Net_in,
                         top_up_net,
                         trial_net_cov_in,
                         vals_pred, 
                         pred_prev,
                         n_mcmc,
                         quant
){
  
  # for the models with the top up of nets there is no trial_net_cov_in column
  # unless the nets are topped up with the trial nets or another type of net there is no top_up_net column
  
  inds <- which(vals_pred$Location_in == Location_in &
                  vals_pred$Net_in == Net_in &
                  vals_pred$top_up_name == top_up_name &
                  vals_pred$Country_in == Country_in &
                  vals_pred$r_model == r_model &
                  vals_pred$int_Net_in == int_Net_in &
                  vals_pred$top_up_net == top_up_net &
                  is.na(vals_pred$trial_net_cov_in) == is.na(trial_net_cov_in))
  
  inds_p <- which(vals_pred$Location_in == Location_in_p &
                  vals_pred$Net_in == "Pyrethroid_only" &
                  vals_pred$top_up_name == top_up_name_p &
                  vals_pred$Country_in == Country_in &
                  vals_pred$r_model == r_model &
                  vals_pred$int_Net_in == "Pyrethroid_only" &
                  vals_pred$top_up_net == "Pyrethroid_only" &
                  is.na(vals_pred$trial_net_cov_in) == is.na(trial_net_cov_in))
  
  if(length(inds) != n_mcmc){return(NA)}else{
    
    prev_mat <- if(Country_in == "Tanzania"){
      sapply(inds, function(i, p){summary_pfpr_0.5_14(p[[i]])}, p = pred_prev)
    } else if(Country_in == "Benin"){
      sapply(inds, function(i, p){summary_pfpr_0_100_all(p[[i]])}, p = pred_prev)
    }
    
    prev_mat_p <- if(Country_in == "Tanzania"){
      sapply(inds_p, function(i, p){summary_pfpr_0.5_14(p[[i]])}, p = pred_prev)
    } else if(Country_in == "Benin"){
      sapply(inds_p, function(i, p){summary_pfpr_0_100_all(p[[i]])}, p = pred_prev)
    }
    
    efficacy <- (prev_mat_p - prev_mat) / prev_mat_p
    
    # returns the mean value from the stochastic simulations
    return(apply(efficacy, 1, quantile, probs = c(quant)))
  }
}

extract_overall_efficacy <- function(Location_in,
                             Net_in,
                             top_up_name,
                             Country_in,
                             r_model,
                             int_Net_in,
                             top_up_net,
                             trial_net_cov_in,
                             vals_pred, 
                             pred_prev,
                             n_mcmc,
                             quant
){
  
  # for the models with the top up of nets there is no trial_net_cov_in column
  # unless the nets are topped up with the trial nets or another type of net there is no top_up_net column
  
  inds <- which(vals_pred$Location_in == Location_in &
                  vals_pred$Net_in == Net_in &
                  vals_pred$top_up_name == top_up_name &
                  vals_pred$Country_in == Country_in &
                  vals_pred$r_model == r_model &
                  vals_pred$int_Net_in == int_Net_in &
                  vals_pred$top_up_net == top_up_net &
                  is.na(vals_pred$trial_net_cov_in) == is.na(trial_net_cov_in))
  
  inds_p <- which(vals_pred$Location_in == Location_in &
                  vals_pred$Net_in == Net_in &
                  vals_pred$top_up_name == top_up_name &
                  vals_pred$Country_in == Country_in &
                  vals_pred$r_model == r_model &
                  vals_pred$int_Net_in == int_Net_in &
                  vals_pred$top_up_net == "none" &
                  vals_pred$trial_net_cov_in == 0)
  
  
  if(length(inds) != n_mcmc){return(NA)}else{
    
    prev_mat <- if(Country_in == "Tanzania"){
      sapply(inds, function(i, p){summary_pfpr_0.5_14(p[[i]])}, p = pred_prev)
    } else if(Country_in == "Benin"){
      sapply(inds, function(i, p){summary_pfpr_0_100_all(p[[i]])}, p = pred_prev)
    }
    
    prev_mat_p <- if(Country_in == "Tanzania"){
      sapply(inds_p, function(i, p){summary_pfpr_0.5_14(p[[i]])}, p = pred_prev)
    } else if(Country_in == "Benin"){
      sapply(inds_p, function(i, p){summary_pfpr_0_100_all(p[[i]])}, p = pred_prev)
    }
    
    efficacy <- (prev_mat_p - prev_mat) / prev_mat_p
    
    # returns the mean value from the stochastic simulations
    return(apply(efficacy, 1, quantile, probs = c(quant)))
  }
}

# point prevalence values at the given times
extract_table <- function(vals_pred,
                          pred_prev,
                          Location_in,
                          Location_in_p = NULL,
                          Net_in,
                          top_up_name,
                          top_up_name_p = NULL,
                          Country_in
){
  
  out <- data.frame(date = pred_prev[[1]]$date) %>% 
    mutate(prev = extract_prev(Location_in = Location_in,
                               Net_in = Net_in, 
                               top_up_name = top_up_name,
                               Country_in = Country_in,
                               r_model = "ll",
                               int_Net_in = Net_in,
                               top_up_net = "Pyrethroid_only",
                               trial_net_cov_in = NA,
                               vals_pred = vals_pred, 
                               pred_prev = pred_prev,
                               n_mcmc = n_mcmc,
                               quant = 0.5),
           
           pred_efficacy = extract_efficacy(Location_in = Location_in,
                                        Location_in_p = Location_in_p,
                                        top_up_name_p = top_up_name_p,
                                        Net_in = Net_in, 
                                        top_up_name = top_up_name,
                                        Country_in = Country_in,
                                        r_model = "ll",
                                        int_Net_in = Net_in,
                                        top_up_net = "Pyrethroid_only",
                                        trial_net_cov_in = NA,
                                        vals_pred = vals_pred, 
                                        pred_prev = pred_prev,
                                        n_mcmc = n_mcmc,
                                        quant = 0.5),
           
           overall_efficacy = extract_overall_efficacy(Location_in = Location_in,
                                                       Net_in = Net_in, 
                                                       top_up_name = top_up_name,
                                                       Country_in = Country_in,
                                                       r_model = "ll",
                                                       int_Net_in = Net_in,
                                                       top_up_net = "Pyrethroid_only",
                                                       trial_net_cov_in = NA,
                                                       vals_pred = vals_pred, 
                                                       pred_prev = pred_prev,
                                                       n_mcmc = n_mcmc,
                                                       quant = 0.5),
           
           tu_prev = extract_prev(Location_in = Location_in,
                                  Net_in = Net_in, 
                                  top_up_name = top_up_name,
                                  Country_in = Country_in,
                                  r_model = "ll",
                                  int_Net_in = Net_in,
                                  top_up_net = Net_in,
                                  trial_net_cov_in = NA,
                                  vals_pred = vals_pred, 
                                  pred_prev = pred_prev,
                                  n_mcmc = n_mcmc,
                                  quant = 0.5),
           
           tu_pred_efficacy = extract_efficacy(Location_in = Location_in,
                                            Location_in_p = Location_in_p,
                                            Net_in = Net_in, 
                                            top_up_name = top_up_name,
                                            top_up_name_p = top_up_name_p,
                                            Country_in = Country_in,
                                            r_model = "ll",
                                            int_Net_in = Net_in,
                                            top_up_net = Net_in,
                                            trial_net_cov_in = NA,
                                            vals_pred = vals_pred, 
                                            pred_prev = pred_prev,
                                            n_mcmc = n_mcmc,
                                            quant = 0.5),
           
           tu_overall_efficacy = extract_overall_efficacy(Location_in = Location_in,
                                                       Net_in = Net_in, 
                                                       top_up_name = top_up_name,
                                                       Country_in = Country_in,
                                                       r_model = "ll",
                                                       int_Net_in = Net_in,
                                                       top_up_net = Net_in,
                                                       trial_net_cov_in = NA,
                                                       vals_pred = vals_pred, 
                                                       pred_prev = pred_prev,
                                                       n_mcmc = n_mcmc,
                                                       quant = 0.5)
    )  %>% mutate(Net = Net_in, Country = Country_in)
  return(out)
}


prev_table <- rbind(
  extract_table(vals_pred = T_vals_pred,
              pred_prev = pred_prev_T,
              Location_in = "Misungwi",
              Net_in = "IG2",
              top_up_name = "top_up_IG2_T",
              Country_in = "Tanzania",
              Location_in_p = "Misungwi",
              top_up_name_p = "top_up_p_only_T"),
              
      extract_table(vals_pred = T_vals_pred,
                       pred_prev = pred_prev_T,
                       Location_in = "Misungwi",
                    Location_in_p = "Misungwi",
                     top_up_name = "top_up_PBO_T",
                       Net_in = "PBO",
                    top_up_name_p = "top_up_p_only_T",
                       Country_in = "Tanzania"),
      extract_table(vals_pred = T_vals_pred,
                              pred_prev = pred_prev_T,
                              Location_in = "Misungwi",
                              Net_in = "Pyrethroid_only",
                              top_up_name = "top_up_p_only_T",
                              Country_in = "Tanzania",
                    Location_in_p = "Misungwi",
                    top_up_name_p = "top_up_p_only_T"),
  
      extract_table(vals_pred = B_vals_pred,
                     pred_prev = pred_prev_B,
                     Location_in = "Zagnanado",
                     Net_in = "IG2",
                     top_up_name = "top_up_IG2_B",
                     Country_in = "Benin",
                    Location_in_p = "Ouinhi",
                    top_up_name_p = "top_up_p_only_B"),
      extract_table(vals_pred = B_vals_pred,
                     pred_prev = pred_prev_B,
                     Location_in = "Ouinhi",
                     Net_in = "Pyrethroid_only",
                     top_up_name = "top_up_p_only_B",
                     Country_in = "Benin",
                    Location_in_p = "Ouinhi",
                    top_up_name_p = "top_up_p_only_B"))



df_out <- left_join(df[, c("Location", "Net", "Country", "Time_months", "Malaria_prevalence", "date")], 
                  prev_table, 
                  by = c("date", "Net", "Country")) %>% subset(Time_months %in% c(6, 12, 18, 24, 30) & Net!= "RG") %>% 
  mutate(
    overall_efficacy = paste0(round(overall_efficacy, digits = 2)*100,"%"),
    tu_overall_efficacy = paste0(round(tu_overall_efficacy, digits = 2)*100,"%"),
         pred_out = paste0(round(prev, digits = 2), " (",round(pred_efficacy, digits = 2)*100,"%)"),
         tu_pred_out = paste0(round(tu_prev, digits = 2), " (",round(tu_pred_efficacy, digits = 2)*100,"%)")) %>% 
  arrange(Country, Net, Time_months)

write.csv(df_out, file = "data/prev_estimates.csv")

# annual mean prevalence values

extract_prev_mean_annual <- function(Location_in,
                                     Location_in_p,
                                     Net_in,
                                     top_up_name,
                                     top_up_name_p,
                                     Country_in,
                                     r_model,
                                     int_Net_in,
                                     top_up_net,
                                     trial_net_cov_in,
                                     vals_pred, 
                                     pred_prev,
                                     n_mcmc,
                                     quant,
                                     s_date,
                                     e_date){
  
  #p <- if(bioassay_uncertainty == "middle"){0.5} else if(bioassay_uncertainty == "lower"){0.025} else{0.975}
  
  inds <- which(vals_pred$Location_in == Location_in &
                  vals_pred$Net_in == Net_in &
                  vals_pred$top_up_name == top_up_name &
                  vals_pred$Country_in == Country_in &
                  vals_pred$r_model == r_model &
                  vals_pred$int_Net_in == int_Net_in &
                  vals_pred$top_up_net == top_up_net &
                  is.na(vals_pred$trial_net_cov_in) == is.na(trial_net_cov_in))
  
  inds_p <- which(vals_pred$Location_in == Location_in_p &
                  vals_pred$Net_in == "Pyrethroid_only" &
                  vals_pred$top_up_name == top_up_name_p &
                  vals_pred$Country_in == Country_in &
                  vals_pred$r_model == r_model &
                  vals_pred$int_Net_in == "Pyrethroid_only" &
                  vals_pred$top_up_net == "Pyrethroid_only" &
                  is.na(vals_pred$trial_net_cov_in) == TRUE)
  
  inds_none <- which(vals_pred$Location_in == Location_in &
                  vals_pred$Net_in == Net_in &
                  vals_pred$top_up_name == top_up_name &
                  vals_pred$Country_in == Country_in &
                  vals_pred$r_model == r_model &
                  vals_pred$int_Net_in == int_Net_in &
                  vals_pred$top_up_net == "none" &
                  vals_pred$trial_net_cov_in == 0)
  
  if(length(inds) != n_mcmc){return(NA)}else{
    
    prev_mat <- if(Country_in == "Tanzania"){sapply(inds, function(i, p){summary_pfpr_0.5_14(p[[i]])}, p = pred_prev)
    } else if(Country_in == "Benin"){sapply(inds, function(i, p){summary_pfpr_0_100_all(p[[i]])}, p = pred_prev)
    }
    
    prev_mat_p <- if(Country_in == "Tanzania"){sapply(inds_p, function(i, p){summary_pfpr_0.5_14(p[[i]])}, p = pred_prev)
    } else if(Country_in == "Benin"){sapply(inds_p, function(i, p){summary_pfpr_0_100_all(p[[i]])}, p = pred_prev)
    }
    
    prev_mat_none <- if(Country_in == "Tanzania"){sapply(inds_none, function(i, p){summary_pfpr_0.5_14(p[[i]])}, p = pred_prev)
    } else if(Country_in == "Benin"){sapply(inds_none, function(i, p){summary_pfpr_0_100_all(p[[i]])}, p = pred_prev)
    }
    
    s_time <- which(pred_prev[[1]]$date == s_date)
    e_time <- which(pred_prev[[1]]$date == e_date)
    
    
    prev_mat <- apply(prev_mat[s_time:e_time,], 2, mean)
    prev_mat_p <- apply(prev_mat_p[s_time:e_time,], 2, mean)
    prev_mat_none <- apply(prev_mat_none[s_time:e_time,], 2, mean)
    
    efficacy <- (prev_mat_p - prev_mat) / prev_mat_p
    overall_efficacy <- (prev_mat_none - prev_mat) / prev_mat_none
    
    return(data.frame(prev = round(quantile(prev_mat, probs = c(quant)), digits = 2),
                      efficacy = paste0(round(quantile(efficacy, probs = c(quant)), digits = 2) * 100, "%"),
                      overall_efficacy = paste0(round(quantile(overall_efficacy, quant), digits = 2) * 100, "%")))
  }
}

extract_prev_mean_annual(Location_in = "Zagnanado",
                         Location_in_p = "Ouinhi",
                         Net_in = "IG2", 
                         top_up_name = "top_up_IG2_B",
                         top_up_name_p = "top_up_p_only_B",
                         Country_in = "Benin",
                         r_model = "ll",
                         int_Net_in = "IG2",
                         top_up_net = "IG2",
                         trial_net_cov_in = NA,
                         vals_pred = B_vals_pred, 
                         pred_prev = pred_prev_B,
                         n_mcmc = n_mcmc,
                         quant = 0.5,
                         s_date = intervention_date_B,
                         e_date = intervention_date_B + 30 * 18)

mean_prev_out <- rbind(data.frame(s_date = rep(c(intervention_date_T, 
                                                 intervention_date_T), 3),
                                  e_date = rep(c(intervention_date_T + 365 * 2,
                                                 intervention_date_T + 365 * 3),3),
                                  year = rep(c("0_to_24", "0_to_36"), 3),
                                  net = sort(rep(c("Pyrethroid_only", "PBO", "IG2"), 2))) %>% mutate(Country = "Tanzania"),
                       
                       data.frame(s_date = rep(c(intervention_date_B, 
                                                 intervention_date_B), 2),
                                  e_date = rep(c(intervention_date_B + 365 * 2,
                                                 intervention_date_B + 365 * 3), 2),
                                  year = rep(c("0_to_24", "0_to_36"), 2),
                                  net = sort(rep(c("Pyrethroid_only", "IG2"), 2))) %>% mutate(Country = "Benin")) %>% 
  rowwise() %>%
  mutate(prev = extract_prev_mean_annual(Net_in = net, 
                                         int_Net_in = net,
                                         Country_in = Country,
                                         vals_pred = if(Country == "Tanzania"){T_vals_pred}else{B_vals_pred}, 
                                         pred_prev = if(Country == "Tanzania"){pred_prev_T}else{pred_prev_B},
                                         top_up_net = NA,
                                         n_rep = 3,
                                         trial_net_cov_in = NULL,
                                         s_date = s_date,
                                         e_date = e_date),
         
         prev_cf = extract_prev_mean_annual(Net_in = net, 
                                            int_Net_in = net,
                                            Country_in = Country,
                                            vals_pred = if(Country == "Tanzania"){T_vals_pred_tu}else{B_vals_pred_tu}, 
                                            pred_prev = if(Country == "Tanzania"){pred_prev_T_tu}else{pred_prev_B_tu},
                                            top_up_net = NA,
                                            n_rep = 3,
                                            trial_net_cov_in = NULL,
                                            s_date = s_date,
                                            e_date = e_date),
         
         prev_none = extract_prev_mean_annual(Net_in = net, 
                                              int_Net_in = net,
                                              Country_in = Country,
                                              vals_pred = if(Country == "Tanzania"){T_vals_pred_tu_none}else{B_vals_pred_tu_none}, 
                                              pred_prev = if(Country == "Tanzania"){pred_prev_T_tu_none}else{pred_prev_B_tu_none},
                                              top_up_net = NA,
                                              n_rep = 3,
                                              trial_net_cov_in = 0,
                                              s_date = s_date,
                                              e_date = e_date),
         
         prev_p_only = extract_prev_mean_annual(Net_in = "Pyrethroid_only", 
                                                int_Net_in = "Pyrethroid_only",
                                                Country_in = Country,
                                                vals_pred = if(Country == "Tanzania"){T_vals_pred}else{B_vals_pred}, 
                                                pred_prev = if(Country == "Tanzania"){pred_prev_T}else{pred_prev_B},
                                                top_up_net = NA,
                                                n_rep = 3,
                                                trial_net_cov_in = NULL,
                                                s_date = s_date,
                                                e_date = e_date)) %>%
  
  mutate(pred_observed = paste0(round(prev, digits = 2), " (", round((prev_p_only - prev)/prev_p_only, digits = 2)*100, "%)"),
         overall_observed = paste0(round((prev_none - prev)/prev_none, digits = 2)*100, "%"),
         
         pred_tu = paste0(round(prev_cf, digits = 2), " (", round((prev_p_only - prev_cf)/prev_p_only, digits = 2)*100, "%)"),
         overall_tu = paste0(round((prev_none - prev_cf)/prev_none, digits = 2)*100, "%"))

write.csv(mean_prev_out, file = "data/mean_prevalence.csv")

# extracting all the different uncertainties

all_combs <- rbind(expand.grid("Country" = c("Tanzania"),
            "r_model" = c("ll","o"),
            "int_Net_in" = c("IG2", "Pyrethroid_only", "PBO"),
            "month" = c(12, 18, 24, 30)
            ) %>% mutate(Net = int_Net_in,
                         Location_in = "Misungwi",
                         top_up_name = case_when(Net == "IG2" ~ "top_up_IG2_T", 
                                                 Net == "PBO" ~ "top_up_PBO_T", 
                                                 Net == "Pyrethroid_only" ~ "top_up_p_only_T"),
                         top_up_name_p = "top_up_p_only_T",
                         Location_in_p = "Misungwi"),
            
            expand.grid("Country" = c("Benin"),
                        "r_model" = c("ll","o"),
                        "int_Net_in" = c("IG2", "Pyrethroid_only"),
                        "month" = c(6, 18, 30)
            ) %>% mutate(Net = int_Net_in,
                         Location_in = ifelse(Net == "IG2", "Zagnanado", "Ouinhi"),
                         top_up_name = ifelse(Net == "IG2", "top_up_IG2_B", "top_up_p_only_B"),
                         top_up_name_p = "top_up_p_only_B",
                         Location_in_p = "Ouinhi")
            ) %>% arrange(Country, int_Net_in, month, r_model)
  
  

all_combs$date <- df[match(interaction(all_combs$Country, all_combs$Net, all_combs$month), interaction(df$Country, df$Net, df$Time_months)),]$date

# uncertainty in the prevalence values
prev_unc <- lapply(seq(1, nrow(all_combs)),
       function(i, all_combs, pred_prev_T, pred_prev_B, T_vals_pred, B_vals_pred){
         out <- data.frame(date = pred_prev_T[[1]]$date) %>% 
           mutate(med = extract_prev(Location_in = all_combs[i, "Location_in"],
                                     Net_in = all_combs[i, "Net"],
                                     top_up_name = all_combs[i, "top_up_name"],
                                     Country_in = all_combs[i, "Country"],
                                     r_model = all_combs[i, "r_model"],
                                     int_Net_in = all_combs[i, "int_Net_in"],
                                     top_up_net = "Pyrethroid_only",
                                     trial_net_cov_in = NA,
                                     vals_pred = if(all_combs[i, "Country"] == "Tanzania"){T_vals_pred} else{B_vals_pred},
                                     pred_prev = if(all_combs[i, "Country"] == "Tanzania"){pred_prev_T} else{pred_prev_B},
                                     n_mcmc = 50,
                                     quant = 0.5),
           lower = extract_prev(Location_in = all_combs[i, "Location_in"],
                                Net_in = all_combs[i, "Net"],
                                top_up_name = all_combs[i, "top_up_name"],
                                Country_in = all_combs[i, "Country"],
                                r_model = all_combs[i, "r_model"],
                                int_Net_in = all_combs[i, "int_Net_in"],
                                top_up_net = "Pyrethroid_only",
                                trial_net_cov_in = NA,
                                vals_pred = if(all_combs[i, "Country"] == "Tanzania"){T_vals_pred} else{B_vals_pred},
                                pred_prev = if(all_combs[i, "Country"] == "Tanzania"){pred_prev_T} else{pred_prev_B},
                                n_mcmc = 50,
                                quant = 0.025),
           upper = extract_prev(Location_in = all_combs[i, "Location_in"],
                        Net_in = all_combs[i, "Net"],
                        top_up_name = all_combs[i, "top_up_name"],
                        Country_in = all_combs[i, "Country"],
                        r_model = all_combs[i, "r_model"],
                        int_Net_in = all_combs[i, "int_Net_in"],
                        top_up_net = "Pyrethroid_only",
                        trial_net_cov_in = NA,
                        vals_pred = if(all_combs[i, "Country"] == "Tanzania"){T_vals_pred} else{B_vals_pred},
                        pred_prev = if(all_combs[i, "Country"] == "Tanzania"){pred_prev_T} else{pred_prev_B},
                        n_mcmc = 50,
                        quant = 0.975))
         ind <- match(all_combs[i, "date"], out$date)
         
         return(paste0(round(out[ind, "med"], digits = 2), " (", round(out[ind, "lower"], digits = 2), "-", round(out[ind, "upper"], digits = 2),")"))
  },
  pred_prev_T = pred_prev_T, pred_prev_B = pred_prev_B, T_vals_pred = T_vals_pred, B_vals_pred = B_vals_pred, all_combs = all_combs)


prev_efficacy_unc <- lapply(seq(1, nrow(all_combs)),
                   function(i, all_combs, pred_prev_T, pred_prev_B, T_vals_pred, B_vals_pred){
                     out <- data.frame(date = pred_prev_T[[1]]$date) %>% 
                       mutate(med = extract_efficacy(Location_in = all_combs[i, "Location_in"],
                                                     Location_in_p = all_combs[i, "Location_in_p"],
                                                 Net_in = all_combs[i, "Net"],
                                                 top_up_name = all_combs[i, "top_up_name"],
                                                 top_up_name_p = all_combs[i, "top_up_name_p"],
                                                 Country_in = all_combs[i, "Country"],
                                                 r_model = all_combs[i, "r_model"],
                                                 int_Net_in = all_combs[i, "int_Net_in"],
                                                 top_up_net = "Pyrethroid_only",
                                                 trial_net_cov_in = NA,
                                                 vals_pred = if(all_combs[i, "Country"] == "Tanzania"){T_vals_pred} else{B_vals_pred},
                                                 pred_prev = if(all_combs[i, "Country"] == "Tanzania"){pred_prev_T} else{pred_prev_B},
                                                 n_mcmc = 50,
                                                 quant = 0.5),
                              lower = extract_efficacy(Location_in = all_combs[i, "Location_in"],
                                                       Location_in_p = all_combs[i, "Location_in_p"],
                                                   Net_in = all_combs[i, "Net"],
                                                   top_up_name = all_combs[i, "top_up_name"],
                                                   top_up_name_p = all_combs[i, "top_up_name_p"],
                                                   Country_in = all_combs[i, "Country"],
                                                   r_model = all_combs[i, "r_model"],
                                                   int_Net_in = all_combs[i, "int_Net_in"],
                                                   top_up_net = "Pyrethroid_only",
                                                   trial_net_cov_in = NA,
                                                   vals_pred = if(all_combs[i, "Country"] == "Tanzania"){T_vals_pred} else{B_vals_pred},
                                                   pred_prev = if(all_combs[i, "Country"] == "Tanzania"){pred_prev_T} else{pred_prev_B},
                                                   n_mcmc = 50,
                                                   quant = 0.025),
                              upper = extract_efficacy(Location_in = all_combs[i, "Location_in"],
                                                       Location_in_p = all_combs[i, "Location_in_p"],
                                                   Net_in = all_combs[i, "Net"],
                                                   top_up_name = all_combs[i, "top_up_name"],
                                                   top_up_name_p = all_combs[i, "top_up_name_p"],
                                                   Country_in = all_combs[i, "Country"],
                                                   r_model = all_combs[i, "r_model"],
                                                   int_Net_in = all_combs[i, "int_Net_in"],
                                                   top_up_net = "Pyrethroid_only",
                                                   trial_net_cov_in = NA,
                                                   vals_pred = if(all_combs[i, "Country"] == "Tanzania"){T_vals_pred} else{B_vals_pred},
                                                   pred_prev = if(all_combs[i, "Country"] == "Tanzania"){pred_prev_T} else{pred_prev_B},
                                                   n_mcmc = 50,
                                                   quant = 0.975))
                     ind <- match(all_combs[i, "date"], out$date)
                     
                     return(paste0(round(out[ind, "med"], digits = 2)*100, "% (", round(out[ind, "lower"], digits = 2)*100, "%-", round(out[ind, "upper"], digits = 2)*100,"%)"))
                   },
                   pred_prev_T = pred_prev_T, pred_prev_B = pred_prev_B, T_vals_pred = T_vals_pred, B_vals_pred = B_vals_pred, all_combs = all_combs)

prev_overall_efficacy <- lapply(seq(1, nrow(all_combs)),
                   function(i, all_combs, pred_prev_T, pred_prev_B, T_vals_pred, B_vals_pred){
                     out <- data.frame(date = pred_prev_T[[1]]$date) %>% 
                       mutate(med = extract_overall_efficacy(Location_in = all_combs[i, "Location_in"],
                                                 Net_in = all_combs[i, "Net"],
                                                 top_up_name = all_combs[i, "top_up_name"],
                                                 Country_in = all_combs[i, "Country"],
                                                 r_model = all_combs[i, "r_model"],
                                                 int_Net_in = all_combs[i, "int_Net_in"],
                                                 top_up_net = "Pyrethroid_only",
                                                 trial_net_cov_in = NA,
                                                 vals_pred = if(all_combs[i, "Country"] == "Tanzania"){T_vals_pred} else{B_vals_pred},
                                                 pred_prev = if(all_combs[i, "Country"] == "Tanzania"){pred_prev_T} else{pred_prev_B},
                                                 n_mcmc = 50,
                                                 quant = 0.5),
                              lower = extract_overall_efficacy(Location_in = all_combs[i, "Location_in"],
                                                   Net_in = all_combs[i, "Net"],
                                                   top_up_name = all_combs[i, "top_up_name"],
                                                   Country_in = all_combs[i, "Country"],
                                                   r_model = all_combs[i, "r_model"],
                                                   int_Net_in = all_combs[i, "int_Net_in"],
                                                   top_up_net = "Pyrethroid_only",
                                                   trial_net_cov_in = NA,
                                                   vals_pred = if(all_combs[i, "Country"] == "Tanzania"){T_vals_pred} else{B_vals_pred},
                                                   pred_prev = if(all_combs[i, "Country"] == "Tanzania"){pred_prev_T} else{pred_prev_B},
                                                   n_mcmc = 50,
                                                   quant = 0.025),
                              upper = extract_overall_efficacy(Location_in = all_combs[i, "Location_in"],
                                                   Net_in = all_combs[i, "Net"],
                                                   top_up_name = all_combs[i, "top_up_name"],
                                                   Country_in = all_combs[i, "Country"],
                                                   r_model = all_combs[i, "r_model"],
                                                   int_Net_in = all_combs[i, "int_Net_in"],
                                                   top_up_net = "Pyrethroid_only",
                                                   trial_net_cov_in = NA,
                                                   vals_pred = if(all_combs[i, "Country"] == "Tanzania"){T_vals_pred} else{B_vals_pred},
                                                   pred_prev = if(all_combs[i, "Country"] == "Tanzania"){pred_prev_T} else{pred_prev_B},
                                                   n_mcmc = 50,
                                                   quant = 0.975))
                     ind <- match(all_combs[i, "date"], out$date)
                     
                     return(paste0(round(out[ind, "med"], digits = 2)*100, "% (", round(out[ind, "lower"], digits = 2)*100, "%-", round(out[ind, "upper"], digits = 2)*100,"%)"))
                   },
                   pred_prev_T = pred_prev_T, pred_prev_B = pred_prev_B, T_vals_pred = T_vals_pred, B_vals_pred = B_vals_pred, all_combs = all_combs)

all_combs <- all_combs %>% mutate(prev_unc = unlist(prev_unc),
                     prev_efficacy_unc = unlist(prev_efficacy_unc),
                     prev_overall_efficacy = unlist(prev_overall_efficacy),
                     r_model = ifelse(r_model == "o", "binomial", "beta-binomial"),
                     int_Net_in = case_when(int_Net_in == "Pyrethroid_only" ~ "Pyrethroid only",
                                            int_Net_in == "IG2" ~ "Pyrethroid-pyrrole",
                                            int_Net_in == "PBO" ~ "Pyrethroid-PBO"))

write.csv(all_combs[,c("Country", "int_Net_in", "month", "r_model", "prev_unc", "prev_efficacy_unc", "prev_overall_efficacy")], file = "data/prev_estimates_uncertainty.csv")

# extracting the prevalence for the MINT comparison
extract_prev_0_5 <- function(Location_in,
                         Net_in,
                         top_up_name,
                         Country_in,
                         r_model,
                         int_Net_in,
                         top_up_net,
                         trial_net_cov_in,
                         vals_pred, 
                         pred_prev,
                         n_mcmc
){
  
  # for the models with the top up of nets there is no trial_net_cov_in column
  # unless the nets are topped up with the trial nets or another type of net there is no top_up_net column
  inds <- which(vals_pred$Location_in == Location_in &
                  vals_pred$Net_in == Net_in &
                  vals_pred$top_up_name == top_up_name &
                  vals_pred$Country_in == Country_in &
                  vals_pred$r_model == r_model &
                  vals_pred$int_Net_in == int_Net_in &
                  vals_pred$top_up_net == top_up_net &
                  is.na(vals_pred$trial_net_cov_in) == is.na(trial_net_cov_in))
  
  if(length(inds) != n_mcmc){return(NA)}else{
    
    prev_mat <- sapply(inds, function(i, p){summary_pfpr_0_5(p[[i]])}, p = pred_prev)
    
    # returns the mean value from the stochastic simulations
    out <- data.frame("date" = pred_prev[[1]]$date,
               "med" = apply(prev_mat, 1, quantile, probs = c(0.5)),
               "lower" = apply(prev_mat, 1, quantile, probs = c(0.025)),
               "upper" = apply(prev_mat, 1, quantile, probs = c(0.975))) %>% 
      mutate(Net = Net_in,
             Country = Country_in,
             Location = Location_in,
             model = ifelse(r_model == "ll", "beta-binomial", "binomial"))
    
    return(out)
  }
}

all_combs_0_5 <- rbind(expand.grid("Country" = c("Tanzania"),
                                   "r_model" = c("ll","o"),
                                   "int_Net_in" = c("IG2", "Pyrethroid_only", "PBO")) %>% 
                         mutate(Net = int_Net_in,
                                Location_in = "Misungwi",
                                top_up_name = case_when(Net == "IG2" ~ "top_up_IG2_T", 
                                     Net == "PBO" ~ "top_up_PBO_T", 
                                     Net == "Pyrethroid_only" ~ "top_up_p_only_T")),
                       expand.grid("Country" = c("Benin"),
            "r_model" = c("ll","o"),
            "int_Net_in" = c("IG2", "Pyrethroid_only")
) %>% mutate(Net = int_Net_in,
             Location_in = ifelse(Net == "IG2", "Zagnanado", "Ouinhi"),
             top_up_name = ifelse(Net == "IG2", "top_up_IG2_B", "top_up_p_only_B")) %>% 
  arrange(Country, int_Net_in, r_model))

MINT_comp_prev <- lapply(seq(1, nrow(all_combs_0_5)), function(i, pred_prev_T, pred_prev_B, T_vals_pred, B_vals_pred, all_combs_0_5){
  return(extract_prev_0_5(Location_in = all_combs_0_5[i, "Location_in"],
               Net_in = all_combs_0_5[i, "Net"],
               top_up_name = all_combs_0_5[i, "top_up_name"],
               Country_in = all_combs_0_5[i, "Country"],
               r_model = all_combs_0_5[i, "r_model"],
               int_Net_in = all_combs_0_5[i, "int_Net_in"],
               top_up_net = "Pyrethroid_only",
               trial_net_cov_in = NA,
               vals_pred = if(all_combs_0_5[i, "Country"] == "Tanzania"){T_vals_pred} else{B_vals_pred},
               pred_prev = if(all_combs_0_5[i, "Country"] == "Tanzania"){pred_prev_T} else{pred_prev_B},
               n_mcmc = 50))}, 
  pred_prev_T = pred_prev_T, pred_prev_B = pred_prev_B, T_vals_pred = T_vals_pred, B_vals_pred = B_vals_pred, all_combs_0_5 = all_combs_0_5)

saveRDS(MINT_comp_prev, file = "data/sim_prev_0_5.rds")
