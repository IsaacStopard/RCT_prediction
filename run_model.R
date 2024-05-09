rm(list = ls())

source(file = "read_data.R")
source(file = "utils/functions.R"); 
source(file = "utils/run_model_functions.R");

# top up values
list2env(readRDS(file = "parameters/top_up_fits.rds"), envir = globalenv())

#####################################
##### running the malaria model #####
#####################################

# calibrating the model to the baseline malaria prevalence

T_vals <- data.frame("Location_in" = rep("Misungwi", 4),
                     "Net_in" = c("IG2", "PBO", "RG", "Pyrethroid_only"),
                     "top_up_name" = c("top_up_IG2_T", "top_up_PBO_T", "top_up_RG_T", "top_up_p_only_T"),
                     "Country_in" = rep("Tanzania", 4))

T_vals <- rbind(T_vals %>% mutate(r_model = "o"), T_vals %>% mutate(r_model = "ll"))

B_vals <- data.frame("Location_in" = c("Cove", "Zagnanado", "Ouinhi"),
                     "Net_in" = c("RG", "IG2", "Pyrethroid_only"),
                     "top_up_name" = c("top_up_RG_B", "top_up_IG2_B", "top_up_p_only_B"),
                     "Country_in" = rep("Benin", 3))

B_vals <- rbind(B_vals %>% mutate(r_model = "o"),
                B_vals %>% mutate(r_model = "ll"))

cl <- makeCluster(5)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("T_site", "int_time_T", "baseline_time_T", "mass_dist_times_T",
                                   "df", "df_species_in", 
                                   "dat_res_pyr_med", "dat_res_pyr_ll_med", "sim_length", 
                                   "top_up_IG2_T", "top_up_PBO_T", 
                                   "top_up_RG_T", "top_up_p_only_T",
                                   "decline_d0", "decline_r0"))

start_EIR_T <- foreach(i=1:nrow(T_vals),
                       .packages = (.packages())
) %dopar% {
  tryCatch({sim_baseline_EIR(
    Location_in = T_vals[i, "Location_in"],
    Net_in = T_vals[i, "Net_in"],
    Country_in = T_vals[i, "Country_in"],
    sites_in = T_site, # foresite package country file
    top_up_name = T_vals[i, "top_up_name"],
    intervention_time = int_time_T,
    baseline_start_time = baseline_time_T,
    df_in = as.data.frame(df),
    prop_species_in = as.data.frame(df_species_in),
    dat_res_pyr_o_in = dat_res_pyr_med,
    dat_res_pyr_ll_in = dat_res_pyr_ll_med,
    r_model = T_vals[i, "r_model"],
    human_population = 10000,
    year = 365,
    sim_years = sim_length/365,
    previous_mass_dist_times = mass_dist_times_T)},
    error = function(cond){
      return(NA)
    })
}
saveRDS(start_EIR_T, 
        file = "results/start_EIR_T_b.rds")

stopCluster(cl)

cl <- makeCluster(5)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("B_site", "int_time_B", "baseline_time_B", "mass_dist_times_B",
                                   "df", "df_species_in", 
                                   "dat_res_pyr_med", "dat_res_pyr_ll_med",
                                   "sim_length", 
                                   "top_up_IG2_B",
                                   "top_up_RG_B", "top_up_p_only_B",
                                   "decline_d0", "decline_r0"))
start_EIR_B <- foreach(i=1:nrow(B_vals),
                       .packages = (.packages())
) %dopar% {
  tryCatch({sim_baseline_EIR(
    Location_in = B_vals[i, "Location_in"],
    Net_in = B_vals[i, "Net_in"],
    Country_in = B_vals[i, "Country_in"],
    sites_in = B_site, # foresite package country file
    top_up_name = B_vals[i, "top_up_name"],
    r_model = B_vals[i, "r_model"],
    intervention_time = int_time_B,
    baseline_start_time = baseline_time_B,
    df_in = as.data.frame(df),
    prop_species_in = as.data.frame(df_species_in),
    dat_res_pyr_o_in = dat_res_pyr_med,
    dat_res_pyr_ll_in = dat_res_pyr_ll_med,
    human_population = 10000,
    year = 365,
    sim_years = sim_length/365,
    previous_mass_dist_times = mass_dist_times_B)},
    error = function(cond){
      return(NA)
    })
}
saveRDS(start_EIR_B, 
        file = "results/start_EIR_B_b.rds")

stopCluster(cl)

start_EIR_T <- readRDS(file = "results/start_EIR_T_b.rds")

start_EIR_B <- readRDS(file = "results/start_EIR_B_b.rds")

##############################
###### running the model #####
##############################

n_mcmc <- 50
mcmc_samples <- sample(seq(1, 1000, 1), n_mcmc, replace = FALSE)

saveRDS(list("n_mcmc" = n_mcmc, "mcmc_samples" = mcmc_samples), file = "results/mcmc_vals.rds")

T_vals <- T_vals %>% mutate(int_Net_in = Net_in, 
                            start_EIR = unlist(start_EIR_T),
                            top_up_net = "Pyrethroid_only",
                            trial_net_cov_in = NA)

T_vals <- rbind(T_vals, 
                subset(T_vals, Net_in != "Pyrethroid_only") %>% mutate(top_up_net = Net_in),
                T_vals %>% mutate(top_up_net = "none",
                                                 trial_net_cov_in = 0))

# T_vals_pred_u <- rbind(T_vals %>% mutate(int_Net_in = "IG2", start_EIR = unlist(start_EIR_T)),
#                       T_vals %>% mutate(int_Net_in = "PBO", start_EIR = unlist(start_EIR_T)),
#                       T_vals %>% mutate(int_Net_in = "Pyrethroid_only", start_EIR = unlist(start_EIR_T)))
# 
# T_vals_pred_b <- rbind(T_vals_pred_u %>% mutate(bioassay_uncertainty = "middle"),
#                         T_vals_pred_u %>% mutate(bioassay_uncertainty = "lower"),
#                         T_vals_pred_u %>% mutate(bioassay_uncertainty = "upper")) %>% as.data.frame()

T_vals_pred <- bind_rows(
  lapply(seq(1, n_mcmc, 1), function(i, T_vals){
    T_vals %>% mutate(index = i)}, 
    T_vals = T_vals)
  )

T_vals_pred <- subset(T_vals_pred, Net_in != "RG") %>% mutate(mcmc_index = mcmc_samples[index])

# Benin
B_vals <- B_vals %>% mutate(int_Net_in = Net_in, 
                            start_EIR = unlist(start_EIR_B),
                            top_up_net = "Pyrethroid_only",
                            trial_net_cov_in = NA)

B_vals <- rbind(B_vals, 
                subset(B_vals, Net_in != "Pyrethroid_only") %>% mutate(top_up_net = Net_in),
                B_vals %>% mutate(top_up_net = "none",
                                  trial_net_cov_in = 0))

B_vals_pred <- bind_rows(
  lapply(seq(1, n_mcmc, 1), function(i, B_vals){
    B_vals %>% mutate(index = i)}, 
    B_vals = B_vals)
)

B_vals_pred <- subset(B_vals_pred, Net_in != "RG") %>% mutate(mcmc_index = mcmc_samples[index])

saveRDS(list("T_vals" = T_vals, "B_vals" = B_vals, "T_vals_pred" = T_vals_pred, "B_vals_pred" = B_vals_pred), file = "results/model_run_combs.rds")

cl <- makeCluster(5)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("T_site", "int_time_T", "baseline_time_T", 
                                   "df", "df_species_in", "sim_length",
                                   "dat_res_pyr",  "dat_res_pp", "dat_res_pbo",
                                   "dat_res_pyr_ll", "dat_res_pp_ll", "dat_res_pbo_ll",
                                   "top_up_IG2_T", "top_up_PBO_T", "mass_dist_times_T",
                                   "top_up_RG_T", "top_up_p_only_T",
                                   "decline_d0", "decline_r0",
                                   "calc_base_mean",
                                   "sim_forward", "T_vals_pred"))

pred_prev_T <- foreach(i=1:nrow(T_vals_pred),
                       .packages = (.packages())
) %dopar% {
  tryCatch({sim_forward(start_EIR = T_vals_pred[i, "start_EIR"],
                        r_model = T_vals_pred[i, "r_model"],
                        
                        dat_res_pyr_o_in = dat_res_pyr,
                        dat_res_pyr_ll_in = dat_res_pyr_ll,
                        
                        dat_res_pp_o_in = dat_res_pp,
                        dat_res_pp_ll_in = dat_res_pp_ll,
                        
                        dat_res_pbo_o_in = dat_res_pbo,
                        dat_res_pbo_ll_in = dat_res_pbo_ll,
                        
                        Location_in = T_vals_pred[i, "Location_in"],
                        Net_in = T_vals_pred[i, "Net_in"],
                        Country_in = T_vals_pred[i, "Country_in"],
                        int_Net_in = T_vals_pred[i, "int_Net_in"],
                        top_up_name = T_vals_pred[i, "top_up_name"],
                        
                        top_up_net = T_vals_pred[i, "top_up_net"],
                        mcmc_index = T_vals_pred[i, "mcmc_index"],
                        
                        sites_in = T_site,
                        intervention_time = int_time_T,
                        baseline_pred_days = baseline_time_T,
                        previous_mass_dist_times = mass_dist_times_T,
                        trial_net_cov_in = T_vals_pred[i, "trial_net_cov_in"])
  },
  error = function(cond){
    return(NA)
  })
}
saveRDS(pred_prev_T, 
        file = "results/pred_prev_T_b_mcmc.rds")

stopCluster(cl)

pred_prev_T <- readRDS(file = "results/pred_prev_T_b_mcmc.rds")

cl <- makeCluster(5)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("B_site", "int_time_B", "baseline_time_B",  "mass_dist_times_B",
                                   "df", "df_species_in",  "sim_length",
                                   "dat_res_pyr", "dat_res_pp", "dat_res_pbo",
                                   "dat_res_pyr_ll", "dat_res_pp_ll", "dat_res_pbo_ll",
                                   "top_up_IG2_B",
                                   "top_up_RG_B", "top_up_p_only_B",
                                   "decline_d0", "decline_r0",
                                   "calc_base_mean", 
                                   "sim_forward", "B_vals_pred"))

pred_prev_B <- foreach(i=1:nrow(B_vals_pred),
                       .packages = (.packages())
) %dopar% {
  tryCatch({sim_forward(start_EIR = B_vals_pred[i, "start_EIR"],
                        r_model = B_vals_pred[i, "r_model"],
                        
                        Location_in = B_vals_pred[i, "Location_in"],
                        Net_in = B_vals_pred[i, "Net_in"],
                        Country_in = B_vals_pred[i, "Country_in"],
                        int_Net_in = B_vals_pred[i, "int_Net_in"],
                        
                        dat_res_pyr_o_in = dat_res_pyr,
                        dat_res_pyr_ll_in = dat_res_pyr_ll,
                        
                        dat_res_pp_o_in = dat_res_pp,
                        dat_res_pp_ll_in = dat_res_pp_ll,
                        
                        dat_res_pbo_o_in = dat_res_pbo,
                        dat_res_pbo_ll_in = dat_res_pbo_ll,
                        
                        top_up_name = B_vals_pred[i, "top_up_name"],
                        top_up_net = B_vals_pred[i, "top_up_net"],
                        
                        mcmc_index = B_vals_pred[i, "mcmc_index"],
                        
                        sites_in = B_site,
                        intervention_time = int_time_B,
                        baseline_pred_days = baseline_time_B,
                        previous_mass_dist_times = mass_dist_times_B,
                        
                        trial_net_cov_in = B_vals_pred[i, "trial_net_cov_in"])
  },
  error = function(cond){
    return(NA)
  })
}
saveRDS(pred_prev_B,
        file = "results/pred_prev_B_b_mcmc.rds")

stopCluster(cl)

# # top up simulations
# T_vals_tu <- T_vals %>% mutate(top_up_net = Net_in)
# 
# # T_vals_tu_s <- subset(T_vals_tu, Net_in != "RG")
# # 
# # T_vals_pred_tu_s <- rbind(T_vals_tu_s %>% mutate(bioassay_uncertainty = "middle"),
# #                        T_vals_tu_s %>% mutate(bioassay_uncertainty = "lower"),
# #                        T_vals_tu_s %>% mutate(bioassay_uncertainty = "upper")) %>% as.data.frame()
# 
# T_vals_pred_tu <- bind_rows(
#   lapply(seq(1, n_mcmc, 1), function(i, T_vals_tu){
#     T_vals_tu %>% mutate(index = i)}, 
#     T_vals_tu = T_vals_tu)
# )
# 
# cl <- makeCluster(5)
# registerDoParallel(cl)
# clusterExport(cl = cl, varlist = c("T_site", "int_time_T", "baseline_time_T", 
#                                    "df", "df_species_in", 
#                                    "dat_res_pyr", "sim_length", "dat_res_pp", "dat_res_pbo",
#                                    "dat_res_pyr_ll", "dat_res_pp_ll", "dat_res_pbo_ll",
#                                    "pyr_binomial_mcmc", "pyr_beta_bin_mcmc", "pp_binomial_mcmc", "pp_beta_bin_mcmc",
#                                    "pbo_binomial_mcmc", "pbo_beta_bin_mcmc",
#                                    "top_up_IG2_T", "top_up_PBO_T", "mass_dist_times_T",
#                                    "top_up_RG_T", "top_up_p_only_T",
#                                    "decline_d0", "decline_r0",
#                                    "calc_base_mean"))
# 
# pred_prev_T_tu <- foreach(i=1:nrow(T_vals_pred_tu),
#                        .packages = (.packages())
# ) %dopar% {
#   tryCatch({sim_forward(start_EIR = T_vals_pred_tu[i, "start_EIR"],
#                         Location_in = T_vals_pred_tu[i, "Location_in"],
#                         Net_in = T_vals_pred_tu[i, "Net_in"],
#                         Country_in = T_vals_pred_tu[i, "Country_in"],
#                         int_Net_in = T_vals_pred_tu[i, "int_Net_in"],
#                         top_up_name = T_vals_pred_tu[i, "top_up_name"],
#                         top_up_net = T_vals_pred_tu[i, "top_up_net"],
#                         r_model = T_vals_pred_tu[i, "r_model"],
#                         sites_in = T_site,
#                         intervention_time = int_time_T,
#                         baseline_pred_days = baseline_time_T,
#                         mcmc_index =  T_vals_pred_tu[i, "index"],
#                         previous_mass_dist_times = mass_dist_times_T)
#   },
#   error = function(cond){
#     return(NA)
#   })
# }
# saveRDS(pred_prev_T_tu, 
#         file = "results/pred_prev_T_tu_b.rds")
# 
# stopCluster(cl)
# 
# pred_prev_T_tu <- readRDS(file = "results/pred_prev_T_tu_b.rds")
# 
# # model with no top up and 0 coverage of trial nets
# 
# 
# cl <- makeCluster(5)
# registerDoParallel(cl)
# clusterExport(cl = cl, varlist = c("T_site", "int_time_T", "baseline_time_T", 
#                                    "df", "df_species_in", 
#                                    "dat_res_pyr", "sim_length", "dat_res_pp", "dat_res_pbo",
#                                    "dat_res_pyr_ll", "dat_res_pp_ll", "dat_res_pbo_ll",
#                                    "pyr_binomial_mcmc", "pyr_beta_bin_mcmc", "pp_binomial_mcmc", "pp_beta_bin_mcmc",
#                                    "pbo_binomial_mcmc", "pbo_beta_bin_mcmc",
#                                    "top_up_IG2_T", "top_up_PBO_T", "mass_dist_times_T",
#                                    "top_up_RG_T", "top_up_p_only_T",
#                                    "decline_d0", "decline_r0",
#                                    "calc_base_mean"))
# 
# pred_prev_T_tu_none <- foreach(i=1:nrow(T_vals_pred_tu_none),
#                           .packages = (.packages())
# ) %dopar% {
#   tryCatch({sim_forward(start_EIR = T_vals_pred_tu_none[i, "start_EIR"],
#                         Location_in = T_vals_pred_tu_none[i, "Location_in"],
#                         Net_in = T_vals_pred_tu_none[i, "Net_in"],
#                         Country_in = T_vals_pred_tu_none[i, "Country_in"],
#                         int_Net_in = T_vals_pred_tu_none[i, "int_Net_in"],
#                         top_up_name = T_vals_pred_tu_none[i, "top_up_name"],
#                         top_up_net = T_vals_pred_tu_none[i, "top_up_net"],
#                         r_model = T_vals_pred_tu_none[i, "r_model"],
#                         sites_in = T_site,
#                         intervention_time = int_time_T,
#                         baseline_pred_days = baseline_time_T,
#                         mcmc_index = T_vals_pred_tu_none[i, "index"],
#                         previous_mass_dist_times = mass_dist_times_T,
#                         trial_net_cov_in = T_vals_pred_tu_none[i, "trial_net_cov_in"])
#   },
#   error = function(cond){
#     return(NA)
#   })
# }
# saveRDS(pred_prev_T_tu_none, 
#         file = "results/pred_prev_T_tu_none_b.rds")
# 
# stopCluster(cl)
# 
# pred_prev_T_tu_none <- readRDS(file = "results/pred_prev_T_tu_none_b.rds")
# 
# 
# pred_prev_B <- readRDS(file = "results/pred_prev_cf_B_b.rds")
# 
# # counterfactual where nets are replaced with old trial nets
# B_vals_tu <- B_vals %>% mutate(int_Net_in = Net_in, 
#                                start_EIR = unlist(start_EIR_B),
#                                top_up_net = Net_in)
# 
# B_vals_tu_s <- subset(B_vals_tu, Net_in != "RG")
# 
# B_vals_pred_tu_s <- rbind(B_vals_tu_s %>% mutate(bioassay_uncertainty = "middle"),
#                           B_vals_tu_s %>% mutate(bioassay_uncertainty = "lower"),
#                           B_vals_tu_s %>% mutate(bioassay_uncertainty = "upper")) %>% as.data.frame()
# 
# B_vals_pred_tu <- bind_rows(
#   lapply(seq(1, 3, 1), function(i, T_vals_pred_b){
#     T_vals_pred_b %>% mutate(rep = i)}, 
#     T_vals_pred_b = B_vals_pred_tu_s)
# )
# 
# cl <- makeCluster(5)
# registerDoParallel(cl)
# clusterExport(cl = cl, varlist = c("B_site", "int_time_B", "baseline_time_B",  "mass_dist_times_B",
#                                    "df", "df_species_in", 
#                                    "dat_res_pyr", "sim_length", "dat_res_pp", "dat_res_pbo",
#                                    "dat_res_pyr_ll", "dat_res_pp_ll", "dat_res_pbo_ll",
#                                    "top_up_IG2_B",
#                                    "top_up_RG_B", "top_up_p_only_B",
#                                    "decline_d0", "decline_r0",
#                                    "calc_base_mean"))
# pred_prev_B_tu <- foreach(i=1:nrow(B_vals_pred_tu),
#                        .packages = (.packages())
# ) %dopar% {
#   tryCatch({sim_forward(start_EIR = B_vals_pred_tu[i, "start_EIR"],
#                         Location_in = B_vals_pred_tu[i, "Location_in"],
#                         Net_in = B_vals_pred_tu[i, "Net_in"],
#                         Country_in = B_vals_pred_tu[i, "Country_in"],
#                         int_Net_in = B_vals_pred_tu[i, "int_Net_in"],
#                         top_up_name = B_vals_pred_tu[i, "top_up_name"],
#                         top_up_net = B_vals_pred_tu[i, "top_up_net"],
#                         r_model = B_vals_pred_tu[i, "r_model"],
#                         sites_in = B_site,
#                         intervention_time = int_time_B,
#                         baseline_pred_days = baseline_time_B,
#                         bioassay_uncertainty = B_vals_pred_tu[i, "bioassay_uncertainty"],
#                         previous_mass_dist_times = mass_dist_times_B)
#   },
#   error = function(cond){
#     return(NA)
#   })
# }
# saveRDS(pred_prev_B_tu, 
#         file = "results/pred_prev_B_tu_b.rds")
# 
# stopCluster(cl)
# 
# pred_prev_B_tu <- readRDS(file = "results/pred_prev_B_tu_b.rds")
# 
# # counterfactual where no trial nets are added
# B_vals_pred_tu_none <- B_vals_pred_tu %>% mutate(top_up_net = "none",
#                                                  trial_net_cov_in = 0)
#                              
# cl <- makeCluster(5)
# registerDoParallel(cl)
# clusterExport(cl = cl, varlist = c("B_site", "int_time_B", "baseline_time_B", "mass_dist_times_B",
#                                    "df", "df_species_in", 
#                                    "dat_res_pyr", "sim_length", "dat_res_pp", "dat_res_pbo",
#                                    "dat_res_pyr_ll", "dat_res_pp_ll", "dat_res_pbo_ll",
#                                    "top_up_IG2_B",
#                                    "top_up_RG_B", "top_up_p_only_B",
#                                    "decline_d0", "decline_r0",
#                                    "calc_base_mean"))
# pred_prev_B_tu_none <- foreach(i=1:nrow(B_vals_pred_tu_none),
#                           .packages = (.packages())
# ) %dopar% {
#   tryCatch({sim_forward(start_EIR = B_vals_pred_tu_none[i, "start_EIR"],
#                         Location_in = B_vals_pred_tu_none[i, "Location_in"],
#                         Net_in = B_vals_pred_tu_none[i, "Net_in"],
#                         Country_in = B_vals_pred_tu_none[i, "Country_in"],
#                         int_Net_in = B_vals_pred_tu_none[i, "int_Net_in"],
#                         top_up_name = B_vals_pred_tu_none[i, "top_up_name"],
#                         top_up_net = B_vals_pred_tu_none[i, "top_up_net"],
#                         r_model = B_vals_pred_tu_none[i, "r_model"],
#                         sites_in = B_site,
#                         intervention_time = int_time_B,
#                         baseline_pred_days = baseline_time_B,
#                         bioassay_uncertainty = B_vals_pred_tu_none[i, "bioassay_uncertainty"],
#                         previous_mass_dist_times = mass_dist_times_B,
#                         trial_net_cov_in = B_vals_pred_tu_none[i, "trial_net_cov_in"])
#   },
#   error = function(cond){
#     return(NA)
#   })
# }
# saveRDS(pred_prev_B_tu_none, 
#         file = "results/pred_prev_B_tu_none_b.rds")
# 
# stopCluster(cl)
# 
# pred_prev_B_tu_none <- readRDS(file = "results/pred_prev_B_tu_none_b.rds")
# 
# ##### saving all the results
# saveRDS(list(
#   "T_vals" = T_vals,
#   "T_vals_pred" = T_vals_pred,
#   "T_vals_pred_tu" = T_vals_pred_tu,
#   "T_vals_pred_tu_none" = T_vals_pred_tu_none,
#   "B_vals" = B_vals,
#   "B_vals_pred" = B_vals_pred,
#   "B_vals_pred_tu" = B_vals_pred_tu,
#   "B_vals_pred_tu_none" = B_vals_pred_tu_none,
#   "pred_prev_T" = pred_prev_T,
#   "pred_prev_T_tu" = pred_prev_T_tu,
#   "pred_prev_T_tu_none" = pred_prev_T_tu_none,
#   "pred_prev_B" = pred_prev_B,
#   "pred_prev_B_tu" = pred_prev_B_tu,
#   "pred_prev_B_tu_none" = pred_prev_B_tu_none), 
# file = "results/model_runs.rds")
