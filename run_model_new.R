rm(list = ls())

# find proportion of people with bednets

# turning on the development mode for the model with LSM

suppressPackageStartupMessages(library(ggplot2)); library(malariasimulation); library(malariaEquilibrium)
library(reshape2); library(tidyverse); library(readxl); library(rstan); library(pracma)
library(foresite); library(doParallel); library(foreach)

source(file = "functions.R"); source(file = "retention_fit_top_up.R")

## Read in seasonality data
sites_csv <- read.csv("parameters/site_parameters.csv",header=TRUE)
# getting the site files for Benin and Tanzania

B_site <- foresite::BEN
T_site <- foresite::TZA

### RCT nets retention and malaria prevalence data
# Tanzania: https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#supplementaryMaterial
# Benin: https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(22)02319-4/fulltext#supplementaryMaterial

df <- read_excel("data/data_RCT.xlsx", sheet = "df") %>% 
  mutate(Bed_net_use_both = Bed_net_use_both/100,
         Bed_net_use_RCT = Bed_net_use_RCT/100,
         Malaria_prevalence = Malaria_prevalence_per_protocol/100)

df_species <- read_excel("data/data_RCT.xlsx", sheet = "site_params")

## read in the parameter efficacy estimates for nets
dat_res_pyr <- read.csv("parameters/pyrethroid_only_nets.csv",header=TRUE) 
dat_res_pbo <- read.csv("parameters/pyrethroid_pbo_nets.csv",header=TRUE) 
dat_res_pp <- read.csv("parameters/pyrethroid_pyrrole_nets.csv",header=TRUE) 

# proportion of bites indoors and in bed from Sherrard-Smith paper (https://www.pnas.org/doi/abs/10.1073/pnas.1820646116)
bite_params <- read_excel("parameters/PNAS.biting.time.xlsx", sheet = "All data for GLMMs")

# calculating the median parameter values by country
bite_params_df <- bite_params %>% group_by(Country_clean, species) %>% summarise(phi_indoor = median(Indoor_phi_mnHuman),
                                                                                 phi_bed = median(Inbed_phi_mnHuman)) %>% 
  mutate(country = ifelse(Country_clean == "BK", "Burkina_Faso", Country_clean))

################################################
##### calculating the top up net coverages #####
################################################

# simulating until 2025
# include a warm up period of ~40 years
sim_length <- 35 * 365

start_date <- as.Date("01/01/1990", format = "%d/%m/%Y")

# Benin times
# March 20th 2020 is the assumed intervention time - when nets are rolled out
# October to November 2019 - baseline time
intervention_date_B <- as.Date("20/03/2020", format = "%d/%m/%Y")
baseline_start_date_B <- as.Date("01/10/2019", format = "%d/%m/%Y")

int_time_B <- difftime(intervention_date_B, start_date)[[1]]
baseline_time_B <- difftime(baseline_start_date_B, start_date)[[1]]

# Tanzania times
#
intervention_date_T <- as.Date("01/01/2020", format = "%d/%m/%Y")
baseline_start_date_T <- as.Date("01/10/2018", format = "%d/%m/%Y")

int_time_T <- difftime(intervention_date_T, start_date)[[1]]
baseline_time_T <- difftime(baseline_start_date_T, start_date)[[1]]

# estimating the observed decline in trial nets
# changing the nets each day
top_up_IG2_B <- get_top_up_multi(net = "IG2",
                                 df = df,
                                 country = "Benin")

top_up_RG_B <- get_top_up_multi(net = "RG",
                                df = df,
                                country = "Benin")

top_up_p_only_B <- get_top_up_multi(net = "Pyrethroid_only",
                                    df = df,
                                    country = "Benin")

top_up_IG2_T <- get_top_up_multi(net = "IG2",
                                 df = df,
                                 country = "Tanzania")

top_up_RG_T <- get_top_up_multi(net = "RG",
                                df = df,
                                country = "Tanzania")

top_up_PBO_T <- get_top_up_multi(net = "PBO",
                                 df = df,
                                 country = "Tanzania")

top_up_p_only_T <- get_top_up_multi(net = "Pyrethroid_only",
                                    df = df,
                                    country = "Tanzania")

top_up_all <- list("top_up_IG2_B" = top_up_IG2_B,
                   "top_up_RG_B" = top_up_RG_B,
                   "top_up_p_only_B" = top_up_p_only_B,
                   "top_up_IG2_T" = top_up_IG2_T,
                   "top_up_RG_T" = top_up_RG_T,
                   "top_up_PBO_T" = top_up_PBO_T,
                   "top_up_p_only_T" = top_up_p_only_T)

saveRDS(top_up_all, file = "parameters/top_up_fits_multinomial.rds")
top_up_all <- readRDS(file = "parameters/top_up_fits_multinomial.rds")

# old methods 
# changing the nets each day
top_up_IG2_B <- get_top_up(net = "IG2",
                           df = df,
                           country = "Benin",
                           tu_diff_time = 1,
                           sim_length = sim_length,
                           int_bed_net_time = int_time_B)

top_up_RG_B <- get_top_up(net = "RG",
                          df = df,
                          country = "Benin",
                          tu_diff_time = 365/2,
                          sim_length = sim_length,
                          int_bed_net_time = int_time_B)

top_up_p_only_B <- get_top_up(net = "Pyrethroid_only",
                              df = df,
                              country = "Benin",
                              tu_diff_time = 365/2,
                              sim_length = sim_length,
                              int_bed_net_time = int_time_B)

top_up_IG2_T <- get_top_up(net = "IG2",
                           df = df,
                           country = "Tanzania",
                           tu_diff_time = 365/2,
                           sim_length = sim_length,
                           int_bed_net_time = int_time_T)

top_up_RG_T <- get_top_up(net = "RG",
                          df = df,
                          country = "Tanzania",
                          tu_diff_time = 365/2,
                          sim_length = sim_length,
                          int_bed_net_time = int_time_T)

top_up_PBO_T <- get_top_up(net = "PBO",
                           df = df,
                           country = "Tanzania",
                           tu_diff_time = 365/2,
                           sim_length = sim_length,
                           int_bed_net_time = int_time_T)

top_up_p_only_T <- get_top_up(net = "Pyrethroid_only",
                              df = df,
                              country = "Tanzania",
                              tu_diff_time = 365/2,
                              sim_length = sim_length,
                              int_bed_net_time = int_time_T)

top_up_all <- list("top_up_IG2_B" = top_up_IG2_B,
                   "top_up_RG_B" = top_up_RG_B,
                   "top_up_p_only_B" = top_up_p_only_B,
                   "top_up_IG2_T" = top_up_IG2_T,
                   "top_up_RG_T" = top_up_RG_T,
                   "top_up_PBO_T" = top_up_PBO_T,
                   "top_up_p_only_T" = top_up_p_only_T)

saveRDS(top_up_all, file = "parameters/top_up_fits.rds")
top_up_all <- readRDS(file = "parameters/top_up_fits.rds")

top_up_IG2_B <- top_up_all$top_up_IG2_B
top_up_RG_B <- top_up_all$top_up_RG_B
top_up_p_only_B <- top_up_all$top_up_p_only_B
top_up_IG2_T <- top_up_all$top_up_IG2_T
top_up_RG_T <- top_up_all$top_up_RG_T
top_up_PBO_T <- top_up_all$top_up_PBO_T
top_up_p_only_T <- top_up_all$top_up_p_only_T

######################################################################
##### calculating the species proportions and bioassay mortality #####
######################################################################

df_species_in <- df_species %>% rowwise() %>% mutate(
  total_n_species = n_gambiae + n_g_or_c + n_coluzzi + n_funestus,
  prop_gambiae = (n_gambiae + n_g_or_c/2) / total_n_species,
  prop_coluzzi = (n_coluzzi + n_g_or_c/2) / total_n_species,
  prop_funestus = n_funestus / total_n_species,
  bioassay_mortality = weighted.mean(x = c(bioassay_gambiae_mortality_y1, bioassay_gambiae_mortality_y2),
                                     w = c(bioassay_gambiae_n_tested_y1, bioassay_gambiae_n_tested_y2), 
                                     na.rm = TRUE)/100,
  pyrethroid_resistance = 1 - bioassay_mortality
)

#####################################
##### running the malaria model #####
#####################################

# more recent seasonality

get_params <- function(Location,
                       Net,
                       Country,
                       sites,
                       int_time,
                       bite_df = as.data.frame(bite_params_df),
                       prop_species_in = as.data.frame(df_species_in)){
  
  
  index <- which(prop_species_in$Location == Location & prop_species_in$Net == Net & prop_species_in$Country == Country)
  
  region <- prop_species_in[index, "Region"] #dide_no <- prop_species_in[index, "dide_code"] ## dide_no calls what is required from site file
  
  site_index <- which(sites$seasonality$name_1 == region)
  site_index_int <- which(sites$interventions$name_1 == region & sites$intervention$urban_rural == "rural" &
                            sites$interventions$year == max(sites$interventions$year))
  
  # setting up the parameter values
  simparams <- get_parameters(
    list(
      human_population = 10000,
      
      prevalence_rendering_min_ages = c(0, 5,  0, 0.5) * 365, ## Prev in 6 months to 14 years measured
      prevalence_rendering_max_ages = c(5,15,100, 14) * 365,
      
      clinical_incidence_rendering_min_ages = c(0, 5,  0, 0.5) * 365, ## All age clin_inc
      clinical_incidence_rendering_max_ages = c(5,15,100, 14) * 365,
      
      severe_incidence_rendering_min_ages = c(0, 5,  0, 0.5) * 365,
      severe_incidence_rendering_max_ages = c(5,15,100, 14) * 365,
      
      model_seasonality = TRUE, ## Seasonality to match study site inputs [sites_13]
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

sim_baseline_EIR_multi <- function(Location_in,
                                   Net_in,
                                   Country_in,
                                   sites_in, # foresite package country file
                                   
                                   top_up_name,
                                   
                                   intervention_time,
                                   baseline_start_time,
                                   
                                   df_in = as.data.frame(df),
                                   prop_species_in = as.data.frame(df_species_in),
                                   dat_res_pyr_in = dat_res_pyr,
                                   
                                   human_population = 10000,
                                   year = 365,
                                   sim_years = sim_length/365
){
  
  top_up <- get(top_up_name)
  
  sim_length <- sim_years * year
  
  simparams <- get_params(Location = Location_in, 
                          Net = Net_in,
                          Country = Country_in,
                          sites = sites_in,
                          int_time = intervention_time)
  
  baseline_index <- which(df_in$Location == Location_in & df_in$Net == Net_in & df_in$Time_months == 0)
  
  baseline_prev <- df_in[baseline_index, "Malaria_prevalence"]
  
  #baseline_prev_days <- intervention_time - 1 #df_in[baseline_index, "Time_months"] * 30 + 365
  
  summary_pfpr <- function(output, 
                           ind = baseline_start_time,
                           min_age = df_in[baseline_index, "Malaria_prev_min_age"]*365,
                           max_age = df_in[baseline_index, "Malaria_prev_max_age"]*365){
    return(mean(output[ind:(ind+30),paste0('n_detect_',min_age,'_',max_age)] / output[ind:(ind+30),paste0('n_',min_age,'_',max_age)]))
  }
  
  index <- which(prop_species_in$Location == Location_in & prop_species_in$Net == Net_in & prop_species_in$Country == Country_in)
  
  bioassay_mortality <- round(prop_species_in[index, "bioassay_mortality"], digits = 2)
  
  bioassay_index <- which(round(dat_res_pyr_in$bioassay_mortality, digits = 2) == bioassay_mortality)
  dn0 <- dat_res_pyr_in[bioassay_index, "dn0_med"]
  rn0 <- dat_res_pyr_in[bioassay_index, "rn0_med"]
  rnm <- 0.24
  gamman <- dat_res_pyr_in[bioassay_index, "gamman_med"]
  
  bednetparams <- simparams
  
  # # initial coverage
  # # last bed net distribution assumed to be 1 year before
  # retention <- mean(rstan::extract(top_up$fit_base, "retention")$retention)
  # rate <- 1 - exp(-1/(retention*365))
  # 
  # # assumed that baseline nets we implemented at time 0 365 days before trial start and give the same mean
  # init_cov <- calc_base_mean(retention = retention,
  #                            mean = top_up$mean_bed_net_use_both,
  #                            min_t = 0,
  #                            max_t = 365)
  # 
  # if(init_cov > 1){
  #   print(paste0("warning: init_cov is greater than 1, mean coverage of baseline year is ",
  #               round(calc_mean_cov(retention = retention, base = 1, 0, 365), 
  #                     digits = 2)
  #               )
  #         )
  #   init_cov <- 1
  # } else{print(paste0("init_cov is less than 1"))}
  
  retention <- mean(rstan::extract(top_up$fit, "retention")$retention) * 365
  
  # need the proportion that you can randomly replace to give this coverage
  # do once every 5 years to give the same mean coverage in the historic use
  
  start_bed_nets <- difftime(as.Date("01/01/2000", format = "%d/%m/%Y"), as.Date("01/01/1990", format = "%d/%m/%Y"))[[1]]
  (int_time_T - start_bed_nets)
  
  mean(subset(sites$interventions, name_1 == prop_species_in[index, "Region"] & urban_rural == "rural" & year < 2020)$itn_use)
  
  
  site_index_int <- which(sites$interventions$name_1 == region & sites$intervention$ & sites$interventions$year == max(sites$interventions$year))
  
  trial_cov <- mean(rstan::extract(top_up$fit, "base[1]")$`base[1]`)
  tot_cov <- mean(rstan::extract(top_up$fit, "base[2]")$`base[2]`) + trial_cov
  base_cov_after  <- tot_cov - trial_cov
  base_cov <- base_cov_after/(1 - trial_cov)
  
  init_cov <- base_cov / exp(-(1 - exp(-1/retention)) * 365)
  
  init_cov <- ifelse(init_cov > 1, 1, init_cov)
  
  # initial bednet event is at time 1
  bednet_events = data.frame(timestep = c(intervention_time - 365),
                             name=c("baseline_nets")
  )
  
  bednetparams_1 <- set_bednets(bednetparams,
                                
                                timesteps = bednet_events$timestep,
                                
                                coverages = init_cov,
                                
                                retention = retention, # assumed
                                
                                # each row needs to show the efficacy parameter across years (and cols are diff mosquito)
                                # gambiae, coluzzi, arabiensis, funestus, coustani
                                # no resistance assumed for coustani
                                
                                dn0 = matrix(c(dn0, dn0, dn0), nrow = 1, ncol = 3),
                                rn = matrix(c(rn0, rn0, rn0), nrow = 1, ncol = 3),
                                rnm = matrix(c(rnm, rnm, rnm), nrow = 1, ncol = 3),
                                gamman = c(gamman * 365)
  )
  
  tol <- ifelse(baseline_prev > 0.7, 0.075, 0.01)
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

T_vals <- data.frame("Location_in" = rep("Misungwi", 4),
                     "Net_in" = c("IG2", "PBO", "RG", "Pyrethroid_only"),
                     "top_up_name" = c("top_up_IG2_T", "top_up_PBO_T", "top_up_RG_T", "top_up_p_only_T"),
                     "Country_in" = rep("Tanzania", 4))

cl <- makeCluster(4)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("T_site", "int_time_T", "baseline_time_T", "df", "df_species_in", 
                                   "dat_res_pyr", "sim_length", 
                                   "top_up_IG2_T", "top_up_PBO_T", 
                                   "top_up_RG_T", "top_up_p_only_T"))
start_EIR_T <- foreach(i=1:nrow(T_vals),
                       .packages = (.packages())
) %dopar% {
  tryCatch({sim_baseline_EIR_multi(
    Location_in = T_vals[i, "Location_in"],
    Net_in = T_vals[i, "Net_in"],
    Country_in = T_vals[i, "Country_in"],
    sites_in = T_site, # foresite package country file
    top_up_name = T_vals[i, "top_up_name"],
    intervention_time = int_time_T,
    baseline_start_time = baseline_time_T,
    df_in = as.data.frame(df),
    prop_species_in = as.data.frame(df_species_in),
    dat_res_pyr_in = dat_res_pyr,
    human_population = 10000,
    year = 365,
    sim_years = sim_length/365)},
    error = function(cond){
      return(NA)
    })
}
saveRDS(start_EIR_T, 
        file = "data/start_EIR_multinomial_T.rds")

stopCluster(cl)
start_EIR_T <- readRDS(file = "data/start_EIR_multinomial_T.rds")


B_vals <- data.frame("Location_in" = c("Cove", "Zagnanado", "Ouinhi"),
                     "Net_in" = c("RG", "IG2", "Pyrethroid_only"),
                     "top_up_name" = c("top_up_RG_B", "top_up_IG2_B", "top_up_p_only_B"),
                     "Country_in" = rep("Benin", 3))

cl <- makeCluster(3)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("B_site", "int_time_B", "baseline_time_B", "df", "df_species_in", 
                                   "dat_res_pyr", "sim_length", 
                                   "top_up_IG2_B",
                                   "top_up_RG_B", "top_up_p_only_B"))
start_EIR_B <- foreach(i=1:nrow(B_vals),
                       .packages = (.packages())
) %dopar% {
  tryCatch({sim_baseline_EIR_multi(
    Location_in = B_vals[i, "Location_in"],
    Net_in = B_vals[i, "Net_in"],
    Country_in = B_vals[i, "Country_in"],
    sites_in = B_site, # foresite package country file
    top_up_name = B_vals[i, "top_up_name"],
    intervention_time = int_time_B,
    baseline_start_time = baseline_time_B,
    df_in = as.data.frame(df),
    prop_species_in = as.data.frame(df_species_in),
    dat_res_pyr_in = dat_res_pyr,
    human_population = 10000,
    year = 365,
    sim_years = sim_length/365)},
    error = function(cond){
      return(NA)
    })
}
saveRDS(start_EIR_B, 
        file = "data/start_EIR_multinomial_B.rds")

stopCluster(cl)

start_EIR_B <- readRDS(file = "data/start_EIR_multinomial_B.rds")

sim_forward_multi <- function(start_EIR,
                              
                              Location_in,
                              Net_in,
                              Country_in,
                              
                              sites_in,
                              
                              int_Net_in,
                              
                              top_up_name,
                              top_up_net,
                              
                              intervention_time,
                              baseline_pred_days,
                              
                              df_in = as.data.frame(df),
                              prop_species_in = as.data.frame(df_species_in),
                              
                              dat_res_pyr_in = dat_res_pyr,
                              dat_res_pp_in = dat_res_pp,
                              dat_res_pbo_in = dat_res_pbo,
                              
                              human_population = 10000,
                              year = 365,
                              sim_years = sim_length/365,
                              top_up_time = 365/2,
                              top_up_cov = 0){
  
  top_up <- get(top_up_name)
  
  # specifying the retention parameter
  # if no top up then let the nets decline
  retention <- mean(rstan::extract(top_up$fit, "retention")$retention) * 365
  
  # need the proportion that you can randomly replace to give this coverage
  trial_cov <- mean(rstan::extract(top_up$fit, "base[1]")$`base[1]`)
  tot_cov <- mean(rstan::extract(top_up$fit, "base[2]")$`base[2]`) + trial_cov
  base_cov_after  <- tot_cov - trial_cov
  base_cov <- base_cov_after/(1 - trial_cov)
  
  init_cov <- base_cov / exp(-(1 - exp(-1/retention)) * 365)
  
  init_cov <- ifelse(init_cov > 1, 1, init_cov)
  
  #
  sim_length <- sim_years * year
  
  simparams <- get_params(Location = Location_in, 
                          Net = Net_in,
                          Country = Country_in,
                          sites = sites_in,
                          int_time = intervention_time)
  
  baseline_index <- which(df_in$Location == Location_in & 
                            df_in$Net == Net_in & 
                            df_in$Country == Country_in &
                            df_in$Time_months == 0)
  baseline_prev <- df_in[baseline_index, "Malaria_prevalence"]
  
  
  index <- which(prop_species_in$Location == Location_in & prop_species_in$Net == Net_in & prop_species_in$Country == Country_in)
  
  bioassay_mortality <- round(prop_species_in[index, "bioassay_mortality"], digits = 2)
  
  bioassay_index <- which(round(dat_res_pyr_in$bioassay_mortality, digits = 2) == bioassay_mortality)
  dn0 <- dat_res_pyr_in[bioassay_index, "dn0_med"]
  rn0 <- dat_res_pyr_in[bioassay_index, "rn0_med"]
  rnm <- 0.24
  gamman <- dat_res_pyr_in[bioassay_index, "gamman_med"]
  
  dat_int <- if(int_Net_in == "IG2"){dat_res_pp_in} else if(int_Net_in == "PBO"){dat_res_pbo_in}else{dat_res_pyr_in}
  
  int_bioassay_index <- which(round(dat_int$bioassay_mortality, digits = 2) == bioassay_mortality)
  dn0_int <- dat_int[int_bioassay_index, "dn0_med"]
  rn0_int <- dat_int[int_bioassay_index, "rn0_med"]
  gamman_int <- dat_int[int_bioassay_index, "gamman_med"]
  
  dat_top_up <- if(top_up_net == "IG2"){dat_res_pp_in}else if(top_up_net == "PBO"){dat_res_pbo_in}else{dat_res_pyr_in}
  
  tu_bioassay_index <- which(round(dat_top_up$bioassay_mortality, digits = 2) == bioassay_mortality)
  dn0_tu <- dat_top_up[tu_bioassay_index, "dn0_med"]
  rn0_tu <- dat_top_up[tu_bioassay_index, "rn0_med"]
  gamman_tu <- dat_top_up[tu_bioassay_index, "gamman_med"]
  
  bednetparams <- simparams
  
  top_up_times <- seq(intervention_time, sim_length, top_up_time)
  n_top_up <- length(top_up_times)
  
  # initial bednet event is at time 1
  bednet_events = data.frame(timestep = c(intervention_time - 365, 
                                          intervention_time, 
                                          top_up_times),
                             name=c("baseline_nets", "intervention_nets", paste0("top_up_nets: ", seq(1, n_top_up)))
  )
  
  bednetparams_1 <- set_bednets(bednetparams,
                                
                                timesteps = bednet_events$timestep,
                                
                                coverages = c(init_cov, init_cov, rep(top_up_cov, n_top_up)),
                                
                                retention = retention, # assumed
                                
                                # each row needs to show the efficacy parameter across years (and cols are diff mosquito)
                                # gambiae, coluzzi, arabiensis, funestus, coustani
                                # no resistance assumed for coustani
                                
                                dn0 = matrix(c(dn0, dn0_int, rep(dn0_tu, n_top_up), dn0, dn0_int, rep(dn0_tu, n_top_up), dn0, dn0_int, rep(dn0_tu, n_top_up)), nrow = n_top_up + 2, ncol = 3),
                                rn = matrix(c(rn0, rn0_int, rep(rn0_tu, n_top_up), rn0, rn0_int, rep(rn0_tu, n_top_up), rn0, rn0_int, rep(rn0_tu, n_top_up)), nrow = n_top_up + 2, ncol = 3),
                                rnm = matrix(c(rnm, rnm, rep(rnm, n_top_up), rnm, rnm, rep(rnm, n_top_up), rnm, rnm, rep(rnm, n_top_up)), nrow = n_top_up + 2, ncol = 3),
                                gamman = c(gamman, gamman_int, rep(gamman_tu, n_top_up)) * 365
  )
  
  correlationsb1 <- get_correlation_parameters(bednetparams_1)
  correlationsb1$inter_round_rho('bednets', 0) # if one nets are given to people who already have nets
  
  bednetparams_1$timesteps <- sim_length 
  
  simparams_eq <- set_equilibrium(bednetparams_1, start_EIR)
  out <- run_simulation(sim_length, simparams_eq, correlations = correlationsb1)
  return(out)
}

T_vals <- data.frame("Location_in" = rep("Misungwi", 3),
                     "Net_in" = c("IG2", "PBO", "Pyrethroid_only"),
                     "top_up_name" = c("top_up_IG2_T", "top_up_PBO_T", "top_up_p_only_T"),
                     "Country_in" = rep("Tanzania", 3),
                     "start_EIR" = c(start_EIR_T[[1]], start_EIR_T[[2]], start_EIR_T[[4]]),
                     "top_up_cov" = rep(0, 3)
)

cl <- makeCluster(4)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("T_site", "int_time_T", "baseline_time_T", "df", "df_species_in", 
                                   "dat_res_pyr", "sim_length", "dat_res_pp", "dat_res_pbo",
                                   "top_up_IG2_T", "top_up_PBO_T", 
                                   "top_up_RG_T", "top_up_p_only_T"))
pred_prev_T <- foreach(i=1:nrow(T_vals),
                       .packages = (.packages())
) %dopar% {
  tryCatch({sim_forward_multi(start_EIR = T_vals[i, "start_EIR"],
                              Location_in = T_vals[i, "Location_in"],
                              Net_in = T_vals[i, "Net_in"],
                              Country_in = T_vals[i, "Country_in"],
                              int_Net_in = T_vals[i, "Net_in"],
                              top_up_name = T_vals[i, "top_up_name"],
                              top_up_net = "none",
                              top_up_time = 365*3,
                              top_up_cov = T_vals[i, "top_up_cov"],
                              sites_in = T_site,
                              intervention_time = int_time_T,
                              baseline_pred_days = baseline_time_T)
  },
  error = function(cond){
    return(NA)
  })
}
saveRDS(pred_prev_T, 
        file = "data/pred_prev_multinomial_T.rds")

stopCluster(cl)

pred_prev_T <- readRDS(file = "data/pred_prev_multinomial_T.rds")

ggplot() +
  geom_line(data = data.frame("t" = pred_prev_T[[1]]$timestep,
                              "prev" = pred_prev_T[[1]]$n_detect_182.5_5110/pred_prev_T[[1]]$n_182.5_5110),
            aes(x = t, y = prev),
            col = "aquamarine") + 
  geom_point(data = subset(df, Location == T_vals[1, "Location_in"] & Net == T_vals[1, "Net_in"]),
             aes(x = Time_months * 30 + baseline_time_T, y = Malaria_prevalence),
             col = "aquamarine") +
  
  geom_line(data = data.frame("t" = pred_prev_T[[2]]$timestep,
                              "prev" = pred_prev_T[[2]]$n_detect_182.5_5110/pred_prev_T[[2]]$n_182.5_5110),
            aes(x = t, y = prev),
            col = "purple") + 
  geom_point(data = subset(df, Location == T_vals[2, "Location_in"] & Net == T_vals[2, "Net_in"]),
             aes(x = Time_months * 30 + baseline_time_T, y = Malaria_prevalence),
             col = "purple") +
  
  geom_line(data = data.frame("t" = pred_prev_T[[3]]$timestep,
                              "prev" = pred_prev_T[[3]]$n_detect_182.5_5110/pred_prev_T[[3]]$n_182.5_5110),
            aes(x = t, y = prev)) + 
  geom_point(data = subset(df, Location == T_vals[3, "Location_in"] & Net == T_vals[3, "Net_in"]),
             aes(x = Time_months * 30 + baseline_time_T, y = Malaria_prevalence),
             col = "grey70", shape = 21, fill = "black") +
  
  coord_cartesian(xlim = c(baseline_time_T - 365*2, 
                           baseline_time_T + 365 * 3)) +
  ylim(c(0, 0.8)) + theme_classic()


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
                             dat_res_pyr_in = dat_res_pyr,
                             
                             human_population = 10000,
                             year = 365,
                             sim_years = sim_length/365
){
  
  top_up <- get(top_up_name)
  
  sim_length <- sim_years * year
  
  simparams <- get_params(Location = Location_in, 
                          Net = Net_in,
                          Country = Country_in,
                          sites = sites_in,
                          int_time = intervention_time)
  
  baseline_index <- which(df_in$Location == Location_in & df_in$Net == Net_in & df_in$Time_months == 0)
  
  baseline_prev <- df_in[baseline_index, "Malaria_prevalence"]
  
  #baseline_prev_days <- intervention_time - 1 #df_in[baseline_index, "Time_months"] * 30 + 365
  
  summary_pfpr <- function(output, 
                           ind = baseline_start_time,
                           min_age = df_in[baseline_index, "Malaria_prev_min_age"]*365,
                           max_age = df_in[baseline_index, "Malaria_prev_max_age"]*365){
    return(mean(output[ind:(ind+30),paste0('n_detect_',min_age,'_',max_age)] / output[ind:(ind+30),paste0('n_',min_age,'_',max_age)]))
  }
  
  index <- which(prop_species_in$Location == Location_in & prop_species_in$Net == Net_in & prop_species_in$Country == Country_in)
  
  bioassay_mortality <- round(prop_species_in[index, "bioassay_mortality"], digits = 2)
  
  bioassay_index <- which(round(dat_res_pyr_in$bioassay_mortality, digits = 2) == bioassay_mortality)
  dn0 <- dat_res_pyr_in[bioassay_index, "dn0_med"]
  rn0 <- dat_res_pyr_in[bioassay_index, "rn0_med"]
  rnm <- 0.24
  gamman <- dat_res_pyr_in[bioassay_index, "gamman_med"]
  
  bednetparams <- simparams
  
  # # initial coverage
  # # last bed net distribution assumed to be 1 year before
  # retention <- mean(rstan::extract(top_up$fit_base, "retention")$retention)
  # rate <- 1 - exp(-1/(retention*365))
  # 
  # # assumed that baseline nets we implemented at time 0 365 days before trial start and give the same mean
  # init_cov <- calc_base_mean(retention = retention,
  #                            mean = top_up$mean_bed_net_use_both,
  #                            min_t = 0,
  #                            max_t = 365)
  # 
  # if(init_cov > 1){
  #   print(paste0("warning: init_cov is greater than 1, mean coverage of baseline year is ",
  #               round(calc_mean_cov(retention = retention, base = 1, 0, 365), 
  #                     digits = 2)
  #               )
  #         )
  #   init_cov <- 1
  # } else{print(paste0("init_cov is less than 1"))}
  
  retention <- 10^10/365
  
  init_cov <- top_up$mean_bed_net_use_both
  
  # initial bednet event is at time 1
  bednet_events = data.frame(timestep = c(1),
                             name=c("baseline_nets")
  )
  
  bednetparams_1 <- set_bednets(bednetparams,
                                
                                timesteps = bednet_events$timestep,
                                
                                coverages = init_cov,
                                
                                retention = retention * year, # assumed
                                
                                # each row needs to show the efficacy parameter across years (and cols are diff mosquito)
                                # gambiae, coluzzi, arabiensis, funestus, coustani
                                # no resistance assumed for coustani
                                
                                dn0 = matrix(c(dn0, dn0, dn0), nrow = 1, ncol = 3),
                                rn = matrix(c(rn0, rn0, dn0), nrow = 1, ncol = 3),
                                rnm = matrix(c(rnm, rnm, rnm), nrow = 1, ncol = 3),
                                gamman = c(gamman * 365)
  )
  
  tol <- ifelse(baseline_prev > 0.7, 0.075, 0.01)
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

start_EIR_IG2_T <- sim_baseline_EIR(Location_in = "Misungwi",
                                    Net_in = "IG2",
                                    Country_in = "Tanzania",
                                    sites_in = T_site, # foresite package country file
                                    
                                    top_up_name = "top_up_IG2_T",
                                    
                                    intervention_time = int_time_T,
                                    baseline_start_time = baseline_time_T)

start_EIR_PBO_T <- sim_baseline_EIR(Location_in = "Misungwi",
                                    Net_in = "PBO",
                                    Country_in = "Tanzania",
                                    sites_in = T_site, # foresite package country file
                                    
                                    top_up_name = "top_up_PBO_T",
                                    
                                    intervention_time = int_time_T,
                                    baseline_start_time = baseline_time_T)

start_EIR_RG_T <- sim_baseline_EIR(Location_in = "Misungwi",
                                   Net_in = "RG",
                                   Country_in = "Tanzania",
                                   sites_in = T_site, # foresite package country file
                                   
                                   top_up_name = "top_up_RG_T",
                                   
                                   intervention_time = int_time_T,
                                   baseline_start_time = baseline_time_T)

start_EIR_p_only_T <- sim_baseline_EIR(Location_in = "Misungwi",
                                       Net_in = "Pyrethroid_only",
                                       Country_in = "Tanzania",
                                       sites_in = T_site, # foresite package country file
                                       
                                       top_up_name = "top_up_p_only_T",
                                       
                                       intervention_time = int_time_T,
                                       baseline_start_time = baseline_time_T)



start_EIR_RG_B <- sim_baseline_EIR(Location_in = "Cove", 
                                   Net_in = "RG", 
                                   Country_in = "Benin",
                                   sites_in = B_site,
                                   top_up_name = "top_up_RG_B",
                                   intervention_time = int_time_B,
                                   baseline_start_time = baseline_time_B)

start_EIR_IG2_B <- sim_baseline_EIR(Location_in = "Zagnanado", 
                                    Net_in = "IG2", 
                                    Country_in = "Benin",
                                    sites_in = B_site,
                                    top_up_name = "top_up_IG2_B",
                                    intervention_time = int_time_B,
                                    baseline_start_time = baseline_time_B)


start_EIR_p_only_B <- sim_baseline_EIR(Location_in = "Ouinhi", 
                                       Net_in = "Pyrethroid_only",
                                       Country_in = "Benin",
                                       sites_in = B_site,
                                       top_up_name = "top_up_p_only_B",
                                       intervention_time = int_time_B,
                                       baseline_start_time = baseline_time_B)


saveRDS(list("start_EIR_IG2_T" = start_EIR_IG2_T, 
             "start_EIR_p_only_T" = start_EIR_p_only_T,
             "start_EIR_RG_T" = start_EIR_RG_T,
             "start_EIR_PBO_T" = start_EIR_PBO_T,
             "start_EIR_IG2_B" = start_EIR_IG2_B,
             "start_EIR_RG_B" = start_EIR_RG_B,
             "start_EIR_p_only_B" = start_EIR_p_only_B), 
        file = "data/start_EIR.rds")

start_EIRs <- readRDS(file = "data/start_EIR_Benin")
start_EIR_IG2 <- start_EIRs$start_EIR_IG2
start_EIR_p_only <- start_EIRs$start_EIR_p_only

sim_forward <- function(start_EIR,
                        
                        Location_in,
                        Net_in,
                        Country_in,
                        
                        int_Net_in,
                        
                        top_up_name,
                        top_up_net,
                        
                        intervention_time,
                        baseline_pred_days,
                        
                        df_in = as.data.frame(df),
                        prop_species_in = as.data.frame(df_species_in),
                        
                        dat_res_pyr_in = dat_res_pyr,
                        dat_res_pp_in = dat_res_pp,
                        dat_res_pbo_in = dat_res_pbo,
                        
                        human_population = 10000,
                        year = 365,
                        sim_years = sim_length/365){
  
  top_up <- get(top_up_name)
  
  # specifying the retention parameter
  # if no top up then let the nets decline
  retention <- if(top_up_net == "none"){mean(rstan::extract(top_up$fit_base, "retention")$retention)}else{10^10/365}
  
  sim_length <- sim_years * year
  
  simparams <- get_params(Location = Location_in, 
                          Net = Net_in,
                          Country = Country_in,
                          sites = sites_in,
                          int_time = intervention_time)
  
  baseline_index <- which(df_in$Location == Location_in & 
                            df_in$Net == Net_in & 
                            df_in$Country == Country_in &
                            df_in$Time_months == 0)
  baseline_prev <- df_in[baseline_index, "Malaria_prevalence"]
  
  
  index <- which(prop_species_in$Location == Location_in & prop_species_in$Net == Net_in & Country == Country_in)
  
  bioassay_mortality <- round(prop_species_in[index, "bioassay_mortality"], digits = 2)
  
  bioassay_index <- which(round(dat_res_pyr_in$bioassay_mortality, digits = 2) == bioassay_mortality)
  dn0 <- dat_res_pyr_in[bioassay_index, "dn0_med"]
  rn0 <- dat_res_pyr_in[bioassay_index, "rn0_med"]
  rnm <- 0.24
  gamman <- dat_res_pyr_in[bioassay_index, "gamman_med"]
  
  dat_int <- if(int_Net_in == "IG2"){dat_res_pp_in} else if(int_Net_in == "PBO"){dat_res_pbo_in}else{dat_res_pyr_in}
  
  int_bioassay_index <- which(round(dat_int$bioassay_mortality, digits = 2) == bioassay_mortality)
  dn0_int <- dat_int[int_bioassay_index, "dn0_med"]
  rn0_int <- dat_int[int_bioassay_index, "rn0_med"]
  gamman_int <- dat_int[int_bioassay_index, "gamman_med"]
  
  dat_top_up <- if(top_up_net == "IG2"){dat_res_pp_in}else if(top_up_net == "PBO"){dat_res_pbo_in}else{dat_res_pyr_in}
  
  tu_bioassay_index <- which(round(dat_top_up$bioassay_mortality, digits = 2) == bioassay_mortality)
  dn0_tu <- dat_top_up[tu_bioassay_index, "dn0_med"]
  rn0_tu <- dat_top_up[tu_bioassay_index, "rn0_med"]
  gamman_tu <- dat_top_up[tu_bioassay_index, "gamman_med"]
  
  bednetparams <- simparams
  
  init_cov <- top_up$mean_bed_net_use_both
  
  if(top_up_net == "none"){
    bednet_events = data.frame(timestep = c(1, covs_df$t + int_time),
                               name=c("baseline_nets", "intervention_nets", paste0("top_up_nets: ", seq(1, nrow(covs_df)-1)))
    )
    
  } else{
    covs_df <- top_up$covs_df
    
    n_top_up <- nrow(covs_df) - 1
    
    # initial bednet event is at time 1
    bednet_events = data.frame(timestep = c(1, covs_df$t + int_time),
                               name=c("baseline_nets", "intervention_nets", paste0("top_up_nets: ", seq(1, nrow(covs_df)-1)))
    )
    
    bednetparams_1 <- set_bednets(bednetparams,
                                  
                                  timesteps = bednet_events$timestep,
                                  
                                  coverages = c(init_cov, init_cov, covs_df[-1, "top_up"]),
                                  
                                  retention = retention, # assumed
                                  
                                  # each row needs to show the efficacy parameter across years (and cols are diff mosquito)
                                  # gambiae, coluzzi, arabiensis, funestus, coustani
                                  # no resistance assumed for coustani
                                  
                                  dn0 = matrix(c(dn0, dn0_int, rep(dn0_tu, n_top_up), dn0, dn0_int, rep(dn0_tu, n_top_up)), nrow = n_top_up + 2, ncol = 2),
                                  rn = matrix(c(rn0, rn0_int, rep(rn0_tu, n_top_up), rn0, rn0_int, rep(rn0_tu, n_top_up)), nrow = n_top_up + 2, ncol = 2),
                                  rnm = matrix(c(rnm, rnm, rep(rnm, n_top_up), rnm, rnm, rep(rnm, n_top_up)), nrow = n_top_up + 2, ncol = 2),
                                  gamman = c(gamman, gamman_int, rep(gamman_tu, n_top_up)) * 365
    )
    
  }
  
  
  
  
  
  
  # # initial coverage
  # # last bed net distribution assumed to be 1 year before
  # retention <- mean(rstan::extract(top_up$fit_base, "retention")$retention)
  # rate <- 1 - exp(-1/(retention*365))
  # 
  # # assumed that baseline nets we implemented at time 0 365 days before trial start and give the same mean
  # init_cov <- calc_base_mean(retention = retention,
  #                            mean = top_up$mean_bed_net_use_both,
  #                            min_t = 0,
  #                            max_t = 365)
  # 
  # if(init_cov > 1){
  #   print(paste0("warning: init_cov is greater than 1, mean coverage of baseline year is ",
  #                round(calc_mean_cov(retention = retention, base = 1, 0, 365), 
  #                      digits = 2)
  #   )
  #   )
  #   init_cov <- 1
  # } else{print(paste0("init_cov is less than 1"))}
  # 
  # n_top_up <- length(top_up$top_up_times)
  # 
  # # initial bednet event is at time 1
  # bednet_events = data.frame(timestep = c(1, 1+ int_bed_net_time, 1+ int_bed_net_time + top_up$top_up_times),
  #                            name=c("baseline_nets", "intervention_nets: ", paste0("top_up_nets", seq(1, n_top_up)))
  # )
  # 
  # 
  # 
  # bednetparams_1 <- set_bednets(bednetparams,
  #                               
  #                               timesteps = bednet_events$timestep,
  #                               
  #                               coverages = c(init_cov, top_up$mean_bed_net_use_both, top_up$top_up_cov),
  #                               
  #                               retention = retention * year, # assumed
  #                               
  #                               # each row needs to show the efficacy parameter across years (and cols are diff mosquito)
  #                               # gambiae, coluzzi, arabiensis, funestus, coustani
  #                               # no resistance assumed for coustani
  #                               
  #                               dn0 = matrix(c(dn0, dn0_int, rep(dn0, n_top_up), dn0, dn0_int, rep(dn0, n_top_up)), nrow = n_top_up + 2, ncol = 2),
  #                               rn = matrix(c(rn0, rn0_int, rep(rn0, n_top_up), rn0, rn0_int, rep(rn0, n_top_up)), nrow = n_top_up + 2, ncol = 2),
  #                               rnm = matrix(c(rnm, rnm, rep(rnm, n_top_up), rnm, rnm, rep(rnm, n_top_up)), nrow = n_top_up + 2, ncol = 2),
  #                               gamman = c(gamman, gamman_int, rep(gamman, n_top_up)) * 365
  # )
  
  correlationsb1 <- get_correlation_parameters(bednetparams_1)
  correlationsb1$inter_round_rho('bednets', 1) # if one nets are given to people who already have nets
  
  bednetparams_1$timesteps <- sim_length 
  
  simparams_eq <- set_equilibrium(bednetparams_1, start_EIR)
  out <- run_simulation(sim_length, simparams_eq, correlations = correlationsb1)
  return(out)
}

out_IG2 <- sim_forward(start_EIR = start_EIR_IG2,
                       Location_in = "Zagnanado",
                       Net_in = "IG2",
                       int_Net_in = "IG2",
                       top_up_name = "top_up_IG2",
                       top_up_net = "pyrethroid_only")

out_IG2_pyr <- sim_forward(start_EIR = start_EIR_IG2,
                           Location_in = "Zagnanado",
                           Net_in = "IG2",
                           int_Net_in = "Pyrethroid_only",
                           top_up_name = "top_up_IG2",
                           top_up_net = "pyrethroid_only")


out_IG2_tu_IG2 <- sim_forward(start_EIR = start_EIR_IG2,
                              Location_in = "Zagnanado",
                              Net_in = "IG2",
                              int_Net_in = "IG2",
                              top_up_name = "top_up_IG2",
                              top_up_net = "IG2")

out_pyrethroid_only <- sim_forward(start_EIR = start_EIR_p_only,
                                   Location_in = "Ouinhi",
                                   Net_in = "Pyrethroid_only",
                                   top_up_name = "top_up_p_only",
                                   top_up_net = "pyrethroid_only")

pred_times <- unique(c(df$Time_months, 36)*30) + int_time
pred_times[1] <- pred_times[1] - 1

pred_indices <- which(out_IG2$timestep %in% pred_times)

round(summary_pfpr_0_100(out_IG2)[pred_indices], digits = 2)

round(summary_pfpr_0_100(out_pyrethroid_only)[pred_indices], digits = 2)

summary_pfpr_0_100_all <- function(output){
  return(output[,'n_detect_0_36500'] / output[,'n_0_36500'])
}


ggplot(data = data.frame("t" = out_pyrethroid_only$timestep,
                         "prev" = summary_pfpr_0_100(out_pyrethroid_only))) +
  geom_line(aes(x = t, y = prev)) +
  geom_point(data = subset(df, Net == "Pyrethroid_only"), 
             aes(x = Time_months * 30 + int_time - 1, y = Malaria_prevalence),
             col = "skyblue", size = 3.5) +
  theme_bw() + 
  xlab("Days") +
  ylab("Malaria prevalence in total population") +
  coord_cartesian(xlim = c(int_time - 365, sim_length))



ggplot(data = data.frame("t" = out_IG2$timestep,
                         "date" = seq(start_date, by = "day", length.out = length(out_IG2$timestep)),
                         "prev" = summary_pfpr_0_100(out_IG2),
                         "prev_top_up_IG2" = summary_pfpr_0_100_all(out_IG2_tu_IG2),
                         "prev_top_up_pyr" = summary_pfpr_0_100_all(out_IG2_pyr))
) +
  geom_line(aes(x = date, y = prev), col = "aquamarine4", alpha = 0.75) +
  geom_line(aes(x = date, y = prev_top_up_IG2), col = "skyblue", alpha = 0.75) +
  geom_line(aes(x = date, y = prev_top_up_pyr), alpha = 0.75) +
  geom_point(data = subset(df, Net == "IG2" & Time_months!=0), 
             aes(x = Time_months * 30 + intervention_date, y = Malaria_prevalence),
             col = "aquamarine4", size = 4.5) +
  geom_segment(x = baseline_start_date, xend = baseline_start_date + 30,
               y = subset(df, Net == "IG2" & Time_months==0)$Malaria_prevalence,
               yend = subset(df, Net == "IG2" & Time_months==0)$Malaria_prevalence,
               col = "aquamarine4", size = 2) +
  geom_vline(xintercept = intervention_date, linetype = 2, size = 1, alpha = 0.5) +
  theme_classic() + 
  xlab("Days") +
  ylab("Malaria prevalence in total population") +
  coord_cartesian(xlim = c(baseline_start_date - 365, baseline_start_date + 365*3))


ggplot(data = data.frame("t" = out_IG2$timestep,
                         "n_use" = out_IG2$n_use_net/10000)) +
  geom_line(aes(x = t, y = n_use)) +
  theme_bw() +
  geom_hline(yintercept = top_up_IG2$mean_bed_net_use_both, linetype = 2) +
  xlab("Days") + ylab("Total net coverage")


out_p_only <- sim_forward(start_EIR = start_EIR_p_only$root,
                          Location_in = "Ouinhi",
                          Net_in = "Pyrethroid_only",
                          top_up_name = "top_up_p_only")

ggplot(data = data.frame("t" = out_p_only$timestep,
                         "prev" = summary_pfpr_0_100(out_p_only))) +
  geom_line(aes(x = t, y = prev)) +
  geom_point(data = subset(df, Net == "Pyrethroid_only"), 
             aes(x = Time_months * 30 + 365 + 1, y = Malaria_prevalence),
             col = "skyblue", size = 3.5) +
  theme_bw() + 
  xlab("Days") +
  ylab("Malaria prevalence in total population")

sum(out_IG2$n_inc_clinical_1825_5475[seq(int_time, int_time + 365)]/out_IG2$n_1825_5475[seq(int_time, int_time + 365)])
out_IG2$clin_