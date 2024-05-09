# functions to help run the malaria model for different trial arms

# more recent seasonality
get_params <- function(Location,
                       Net,
                       Country,
                       sites,
                       int_time,
                       bite_df = as.data.frame(bite_params_df), #
                       prop_species_in = as.data.frame(df_species_in)){ #
  
  
  index <- which(prop_species_in$Location == Location & prop_species_in$Net == Net & prop_species_in$Country == Country)
  
  region <- prop_species_in[index, "Region"] #dide_no <- prop_species_in[index, "dide_code"] ## dide_no calls what is required from site file
  
  site_index <- which(sites$seasonality$name_1 == region)
  site_index_int <- which(sites$interventions$name_1 == region & sites$intervention$urban_rural == "rural" &
                            sites$interventions$year == max(sites$interventions$year))
  
  # setting up the parameter values
  simparams <- get_parameters(
    list(
      human_population = 10000,
      
      prevalence_rendering_min_ages = c(0, 5,  0, 0.5, 0.5) * 365, ## Prev in 6 months to 14 years measured
      prevalence_rendering_max_ages = c(5,15,100, 14, 10) * 365,
      
      clinical_incidence_rendering_min_ages = c(0, 5,  0, 0.5, 0.5) * 365, ## All age clin_inc
      clinical_incidence_rendering_max_ages = c(5,15,100, 14, 10) * 365,
      
      severe_incidence_rendering_min_ages = c(0, 5,  0, 0.5, 0.5) * 365,
      severe_incidence_rendering_max_ages = c(5,15,100, 14, 10) * 365,
      
      #model_seasonality = TRUE, ## Seasonality to match study site inputs [sites_13]
      model_seasonality = TRUE,
      g0 = sites$seasonality$g0[site_index],
      g = c(sites$seasonality$g1[site_index], sites$seasonality$g2[site_index], sites$seasonality$g3[site_index]),
      h = c(sites$seasonality$h1[site_index], sites$seasonality$h2[site_index], sites$seasonality$h3[site_index]),
      
      individual_mosquitoes = FALSE, ## Update next
      
      bednets = TRUE,
      
      smc = FALSE
    )
  )
  
  # set treatment
  simparams <- set_drugs(simparams, list(AL_params, SP_AQ_params, DHA_PQP_params))
  simparams <- set_clinical_treatment(simparams, 
                                      drug=1,
                                      time=c(1),
                                      coverage=c(sites$interventions[site_index_int,"prop_act"] * sites$interventions[site_index_int,"tx_cov"])) # baseline coverage from foresite package
  simparams <- set_clinical_treatment(simparams, 
                                      drug=2,
                                      time=c(1),
                                      coverage=c(sites$interventions[site_index_int,"tx_cov"] - (sites$interventions[site_index_int,"prop_act"] * sites$interventions[site_index_int,"tx_cov"])))
  simparams <- set_clinical_treatment(simparams, 
                                      drug=3,
                                      time=c(1),
                                      coverage=0)
  
  bite_params_in = subset(bite_df, country == prop_species_in[index, "Country"])
  
  gam_params <- malariasimulation::gamb_params # gambiae
  coluzzi_params <- malariasimulation::gamb_params
  fun_params <- malariasimulation::fun_params
  
  gam_params['species'] <- "gambiae"
  coluzzi_params['species'] <- "coluzzi"
  fun_params['species'] <- "funestus"
  
  gam_params[["phi_bednets"]] <- bite_params_in[bite_params_in$species == "An_gambiae_sl","phi_bed"][[1]]
  gam_params[["phi_indoors"]] <- bite_params_in[bite_params_in$species == "An_gambiae_sl","phi_indoor"][[1]]
  
  # assuming the coluzzi has the same parameters as gambiae
  coluzzi_params[["phi_bednets"]] <- bite_params_in[bite_params_in$species == "An_gambiae_sl","phi_bed"][[1]]
  coluzzi_params[["phi_indoors"]] <- bite_params_in[bite_params_in$species == "An_gambiae_sl","phi_indoor"][[1]]
  
  fun_params[["phi_bednets"]] <- bite_params_in[bite_params_in$species == "An_funestus","phi_bed"][[1]]
  fun_params[["phi_indoors"]] <- bite_params_in[bite_params_in$species == "An_funestus","phi_indoor"][[1]]
  
  simparams <- set_species(simparams,
                           list(gam_params, 
                                coluzzi_params,
                                fun_params),
                           c(prop_species_in[index, "prop_gambiae"], 
                             prop_species_in[index, "prop_coluzzi"],
                             prop_species_in[index, "prop_funestus"]))
  
  return(simparams)
}

# original methods
sim_baseline_EIR <- function(Location_in,
                             Net_in,
                             Country_in,
                             sites_in, # foresite package country file
                             
                             top_up_name,
                             
                             intervention_time,
                             baseline_start_time,
                             
                             df_in = as.data.frame(df),
                             prop_species_in = as.data.frame(df_species_in),
                             dat_res_pyr_o_in,
                             dat_res_pyr_ll_in,
                             
                             r_model,
                             
                             human_population = 10000,
                             year = 365,
                             sim_years = sim_length/365,
                             previous_mass_dist_times
){
  
  dat_res_pyr_in <- if(r_model == "o"){dat_res_pyr_o_in}else{dat_res_pyr_ll_in}
  
  top_up <- get(top_up_name)
  
  sim_length <- sim_years * year
  
  simparams <- get_params(Location = Location_in, 
                          Net = Net_in,
                          Country = Country_in,
                          sites = sites_in,
                          int_time = intervention_time)
  
  baseline_index <- which(df_in$Location == Location_in & df_in$Net == Net_in & is.na(df_in$Time_months) == 1)
  
  baseline_prev <- df_in[baseline_index, "Malaria_prevalence"]
  
  #baseline_prev_days <- intervention_time - 1 #df_in[baseline_index, "Time_months"] * 30 + 365
  
  # done for the 30 days
  summary_pfpr <- function(output, 
                           ind = baseline_start_time,
                           min_age = df_in[baseline_index, "Malaria_prev_min_age"]*365,
                           max_age = df_in[baseline_index, "Malaria_prev_max_age"]*365){
    return(mean(output[(ind-15):(ind+15),paste0('n_detect_',min_age,'_',max_age)] / output[(ind-15):(ind+15),paste0('n_',min_age,'_',max_age)]))
  }
  
  index <- which(prop_species_in$Location == Location_in & prop_species_in$Net == Net_in & prop_species_in$Country == Country_in)
  
  bioassay_mortality <- round(prop_species_in[index, "bioassay_mortality"], digits = 2)
  
  bioassay_index <- which(round(dat_res_pyr_in$bioassay_mortality, digits = 2) == bioassay_mortality)
  dn0 <- dat_res_pyr_in[bioassay_index, "dn0_med"]
  rn0 <- dat_res_pyr_in[bioassay_index, "rn0_med"]
  rnm <- 0.24
  gamman <- dat_res_pyr_in[bioassay_index, "gamman_med"]
  
  bednetparams <- simparams
  
  retention <- 10^10
  
  init_cov <- min(df_in[baseline_index, "Bed_net_use_both"], top_up$mean_bed_net_use_both) # # changed to baseline coverage
  
  # initial bednet event is at time 1
  n_mass_dist <- length(previous_mass_dist_times)
  
  bednet_events = data.frame(timestep = c(1, previous_mass_dist_times),
                             name=c("baseline_nets", paste0("mass distribution: ", seq(1, n_mass_dist)))
  )
  
  bednetparams_1 <- set_bednets(bednetparams,
                                
                                timesteps = bednet_events$timestep,
                                
                                coverages = rep(init_cov, (1 + n_mass_dist)),
                                
                                retention = retention, # assumed
                                
                                dn0 = matrix(c(dn0, dn0, dn0,
                                               rep(dn0, n_mass_dist * 3)), nrow = (1 + n_mass_dist), ncol = 3),
                                rn = matrix(c(rn0, rn0, rn0,
                                              rep(rn0, n_mass_dist * 3)), nrow = (1 + n_mass_dist), ncol = 3),
                                rnm = matrix(c(rnm, rnm, rnm,
                                               rep(rnm, n_mass_dist * 3)), nrow = (1 + n_mass_dist), ncol = 3),
                                gamman = rep(gamman/log(2) * 365,  (1 + n_mass_dist)) # divided by log(2) because gamman is treated as the mean in malariasim
  )
  
  tol <- ifelse(baseline_prev > 0.7, 0.075, 0.0075)
  set.seed(12345)
  bednetparams_1$timesteps <- sim_length 
  out <- cali::calibrate(parameters = bednetparams_1,
                         target = baseline_prev,
                         #target_tt = baseline_prev_days,
                         summary_function = summary_pfpr,
                         tolerance = tol,
                         weights = 1,
                         low = 0.1, 
                         high = 1500)
  #rm(list = c(top_up, bednetparams_1, simparams, baseline_index, baseline_prev, index))
  return(out)
}

sim_forward <- function(start_EIR,
                        r_model,
                        
                        Location_in,
                        Net_in,
                        Country_in,
                        
                        int_Net_in,
                        
                        top_up_name,
                        
                        top_up_net,
                        mcmc_index,
                        sites_in,
                        intervention_time,
                        baseline_pred_days,
                        
                        df_in = as.data.frame(df),
                        prop_species_in = as.data.frame(df_species_in),
                        
                        dat_res_pyr_o_in,
                        dat_res_pyr_ll_in,
                        
                        dat_res_pp_o_in,
                        dat_res_pp_ll_in,
                        
                        dat_res_pbo_o_in,
                        dat_res_pbo_ll_in,
                        
                        human_population = 10000,
                        year = 365,
                        sim_years = sim_length/365,
                        
                        previous_mass_dist_times,
                        trial_net_cov_in){
  
  print(paste0(Location_in, ": ", int_Net_in))
  
  top_up <- get(top_up_name)
  
  # specifying the retention parameter
  # if no top up then let the nets decline
  retention <- if(top_up_net == "none"){
    median(rstan::extract(top_up$fit_base, "retention")$retention) * 365
  }else{10^10}
  
  sim_length <- sim_years * year
  
  simparams <- get_params(Location = Location_in, 
                          Net = Net_in,
                          Country = Country_in,
                          sites = sites_in,
                          int_time = intervention_time)
  
  baseline_index <- which(df_in$Location == Location_in & 
                            df_in$Net == Net_in & 
                            df_in$Country == Country_in &
                            is.na(df_in$Time_months) == 1)
  baseline_prev <- df_in[baseline_index, "Malaria_prevalence"]
  
  
  index <- which(prop_species_in$Location == Location_in & prop_species_in$Net == Net_in & prop_species_in$Country == Country_in)
  
  bioassay_mortality <- round(prop_species_in[index, "bioassay_mortality"], digits = 2)
  
  # if bioassay mortality is lower then the lower death probability estimate is used else the higher value if used
  
  dat_res_pyr_in <- if(r_model == "o"){dat_res_pyr_o_in}else{dat_res_pyr_ll_in}
  dat_res_pp_in <- if(r_model == "o"){dat_res_pp_o_in}else{dat_res_pp_ll_in}
  dat_res_pbo_in <- if(r_model == "o"){dat_res_pbo_o_in}else{dat_res_pbo_ll_in}
  
  bioassay_index <- which(round(dat_res_pyr_in$bioassay_mortality, digits = 2) == bioassay_mortality & dat_res_pyr_in$index == mcmc_index)
  
  dn0 <- dat_res_pyr_in[bioassay_index, "dn0"]
  rn0 <- dat_res_pyr_in[bioassay_index, "rn0"]
  rnm <- 0.24
  gamman <- dat_res_pyr_in[bioassay_index, "gamman"] * 365
  
  dat_int <- if(int_Net_in == "IG2"){dat_res_pp_in} else if(int_Net_in == "PBO"){dat_res_pbo_in} else if(int_Net_in == "PPF"){dat_res_ppf_in} else if(int_Net_in == "Pyrethroid_only"){dat_res_pyr_in}
  
  int_bioassay_index <- which(round(dat_int$bioassay_mortality, digits = 2) == bioassay_mortality & dat_int$index == mcmc_index)
  
  dn0_int <- dat_int[int_bioassay_index, "dn0"]
  rn0_int <- dat_int[int_bioassay_index, "rn0"]
  gamman_int <- dat_int[int_bioassay_index, "gamman"] * 365
  
  dat_top_up <- if(top_up_net == "IG2"){dat_res_pp_in}else if(top_up_net == "PBO"){dat_res_pbo_in} else if(top_up_net == "PPF"){dat_res_ppf_in} else if(top_up_net == "Pyrethroid_only"){dat_res_pyr_in}
  
  if(top_up_net != "none"){
    tu_bioassay_index <- which(round(dat_top_up$bioassay_mortality, digits = 2) == bioassay_mortality & dat_top_up$index == mcmc_index)
    dn0_tu <- dat_top_up[tu_bioassay_index, "dn0"]
    rn0_tu <- dat_top_up[tu_bioassay_index, "rn0"]
    gamman_tu <- dat_top_up[tu_bioassay_index, "gamman"] * 365
    
  }
  
  bednetparams <- simparams
  
  init_cov <- min(df_in[baseline_index, "Bed_net_use_both"], top_up$mean_bed_net_use_both) # 
  covs_df <- top_up$covs_df
  
  n_mass_dist <- length(previous_mass_dist_times)
  
  # if no replacement of lost nets
  # assumes the same constant coverage prior to the net distribution
  
  if(top_up_net == "none"){
    
    trial_net_cov <- if(is.na(trial_net_cov_in) == FALSE){trial_net_cov_in}else{top_up$mean_bed_net_use_both} # default is to use the mean bed net use value
    
    init_time_diff <- 30
    
    # one extra previous mass distribution
    init_times <- seq(1, (previous_mass_dist_times - init_time_diff), init_time_diff)
    dn0_seq <- decline_d0(t = (init_times - 1), dn0 = dn0, gamma_n = gamman) # - 1 because the initial coverage is at time 1
    rn0_seq <- decline_r0(t = (init_times - 1), rn0 = rn0, rnm = rnm, gamma_n = gamman) # - 1 because the initial coverage is at time 1
    
    placeholder_times <- seq(previous_mass_dist_times, intervention_time - init_time_diff, init_time_diff)
    place_dn0_seq <- decline_d0(t = (placeholder_times - previous_mass_dist_times), dn0 = dn0, gamma_n = gamman)
    place_rn0_seq <- decline_r0(t = (placeholder_times - previous_mass_dist_times), rn0 = rn0, rnm = rnm, gamma_n = gamman)
    
    init_times <- c(init_times, placeholder_times)
    dn0_seq <- c(dn0_seq, place_dn0_seq)
    rn0_seq <- c(rn0_seq, place_rn0_seq)
    
    #init_times <- round(seq(1, intervention_time, 30), digits = 0)
    n_init_times <- length(init_times)
    
    # calculating the initial coverage required to give the mean init cov values
    mean_init_cov <- calc_base_mean(retention = retention, 
                                    mean = init_cov,
                                    min_t = 0,
                                    max_t = init_time_diff)
    
    bednet_events = data.frame(timestep = c(init_times, intervention_time),
                               name=c(paste0(rep("baseline_nets", n_init_times),
                                             ": ", init_times), 
                                      "intervention_nets"))
    
    bednetparams_1 <- set_bednets(bednetparams,
                                  
                                  timesteps = bednet_events$timestep,
                                  
                                  coverages = c(rep(mean_init_cov, n_init_times), trial_net_cov),
                                  
                                  retention = retention, # assumed
                                  
                                  # each row needs to show the efficacy parameter across years (and cols are diff mosquito)
                                  # gambiae, coluzzi, arabiensis, funestus, coustani
                                  # no resistance assumed for coustani
                                  
                                  dn0 = matrix(c(dn0_seq, dn0_int, 
                                                 dn0_seq, dn0_int, 
                                                 dn0_seq, dn0_int), nrow = n_init_times + 1, ncol = 3),
                                  rn = matrix(c(rn0_seq, rn0_int, 
                                                rn0_seq, rn0_int, 
                                                rn0_seq, rn0_int), nrow = n_init_times + 1, ncol = 3),
                                  rnm = matrix(c(rep(rnm, n_init_times), rnm,
                                                 rep(rnm, n_init_times), rnm, 
                                                 rep(rnm, n_init_times), rnm), nrow = n_init_times + 1, ncol = 3),
                                  gamman = c(rep(gamman/log(2), n_init_times), gamman_int/log(2))
    )
    
  } else{
    
    covs_df <- covs_df[-nrow(covs_df),]
    n_top_up <- nrow(covs_df)
    t_times <- covs_df$t + intervention_time
    p_times <- t_times[-1] - 1
    n_p_times <- length(p_times)
    n_t_times <- length(t_times)
    int_times <- na.omit(c(rbind(t_times, c(p_times, NA))))
    
    p_covs <- rep(top_up$mean_bed_net_use_both, n_p_times)
    
    net_names <- na.omit(c(rbind(paste0("Trial nets: ", seq(1, n_t_times)), c(paste0("Replacement nets: ", seq(1, n_p_times)), NA)))) %>% as.vector()
    
    # pyrethroid only nets that are replacing are assumed to be the same age as the trial nets
    ds_int <- decline_d0(covs_df$t, dn0 = dn0_int, gamma_n = gamman_int)
    ds_tu <- decline_d0(p_times - intervention_time, dn0 = dn0_tu, gamma_n = gamman_tu)
    
    rs_int <- decline_r0(covs_df$t, rn0 = rn0_int, gamma_n = gamman_int, rnm = rnm)
    rs_tu <- decline_r0(p_times - intervention_time, rn0 = rn0_tu, gamma_n = gamman_tu, rnm = rnm)
    
    dn0_vals <- na.omit(c(rbind(ds_int, c(ds_tu, NA))))
    rn0_vals <- na.omit(c(rbind(rs_int, c(rs_tu, NA))))
    
    gamman_in <- na.omit(c(rbind(rep(gamman_int, n_t_times), c(rep(gamman_tu, n_p_times), NA)))) %>% as.vector()
    
    top_up_name <- "top_up" # assumes no uncertainty in the net retention times
    
    covs_in <- na.omit(c(rbind(covs_df[,top_up_name], c(p_covs, NA)))) %>% as.vector()
    
    # initial bednet event is at time 1
    
    bednet_events = data.frame(timestep = c(1, previous_mass_dist_times, int_times),
                               name=c("baseline_nets", paste0("mass distribution: ", seq(1, n_mass_dist)), net_names)
    )
    
    dn0_mass_dist <- rep(dn0, n_mass_dist)
    rn0_mass_dist <- rep(rn0, n_mass_dist)
    rnm_mass_dist <- rep(rnm, n_mass_dist)
    
    bednetparams_1 <- set_bednets(bednetparams,
                                  
                                  timesteps = bednet_events$timestep,
                                  
                                  coverages = c(rep(init_cov, (1 + n_mass_dist)), covs_in), #c(init_cov, covs_df[1, init_cov_name], covs_df[-1, top_up_name]), # initial coverage is the mean coverage between first two
                                  
                                  retention = retention, # assumed
                                  
                                  # each row needs to show the efficacy parameter across years (and cols are diff mosquito)
                                  # gambiae, coluzzi, arabiensis, funestus, coustani
                                  # no resistance assumed for coustani
                                  
                                  dn0 = matrix(c(dn0, dn0_mass_dist, dn0_vals,
                                                 dn0, dn0_mass_dist, dn0_vals,
                                                 dn0, dn0_mass_dist, dn0_vals), nrow = length(int_times) + n_mass_dist + 1, ncol = 3),
                                  
                                  rn = matrix(c(rn0, rn0_mass_dist, rn0_vals,
                                                rn0, rn0_mass_dist, rn0_vals,
                                                rn0, rn0_mass_dist, rn0_vals), nrow = length(int_times) + n_mass_dist + 1, ncol = 3),
                                  
                                  rnm = matrix(c(rnm, rnm_mass_dist, rep(rnm, length(int_times)),
                                                 rnm, rnm_mass_dist, rep(rnm, length(int_times)),
                                                 rnm, rnm_mass_dist, rep(rnm, length(int_times))), nrow = length(int_times) + n_mass_dist + 1, ncol = 3),
                                  
                                  gamman = c(gamman/log(2), rep(gamman/log(2), n_mass_dist), gamman_in/log(2)))
    
  }
  
  correlationsb1 <- get_correlation_parameters(bednetparams_1)
  correlationsb1$inter_round_rho('bednets', 1) # if one nets are given to people who already have nets
  
  bednetparams_1$timesteps <- sim_length 
  
  simparams_eq <- set_equilibrium(parameters = bednetparams_1, init_EIR = start_EIR)
  out <- run_simulation(sim_length, simparams_eq, correlations = correlationsb1)
  return(out)
}

# functions to extract the results
summary_pfpr_0_100_all <- function(output){
  return(output[,'n_detect_0_36500'] / output[,'n_0_36500'])
}

summary_pfpr_0.5_14 <- function(output){
  return(output[,'n_detect_182.5_5110'] / output[,'n_182.5_5110'])
}

summary_pfpr_0_5 <- function(output){
  return(output[,'n_detect_0_1825'] / output[,'n_0_1825'])
}

