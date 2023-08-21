rm(list = ls())

# find proportion of people with bednets

# turning on the development mode for the model with LSM

suppressPackageStartupMessages(library(ggplot2)); library(malariasimulation); library(malariaEquilibrium)
library(reshape2); library(tidyverse); library(readxl); library(rstan); library(pracma)
library(foresite); library(doParallel); library(foreach);
library(patchwork); library(pammtools); library(ggpattern);
library(cowplot)

source(file = "functions.R"); source(file = "retention_fit_top_up.R")

### functions

## Read in seasonality data
# getting the site files for Benin and Tanzania

B_site <- foresite::BEN
T_site <- foresite::TZA

### RCT nets retention and malaria prevalence data
# Tanzania: https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#supplementaryMaterial
# Benin: https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(22)02319-4/fulltext#supplementaryMaterial

df <- read_excel("data/data_RCT.xlsx", sheet = "df") %>% 
  mutate(Bed_net_use_both = Bed_net_use_both/100,
         Bed_net_use_RCT = Bed_net_use_RCT/100,
         Malaria_prevalence = Malaria_prevalence_per_protocol/100,
         Malaria_prevalence_l = Malaria_prevalence - 1.96*sqrt(Malaria_prevalence * (1 - Malaria_prevalence) / Malaria_prev_n_tested),
         Malaria_prevalence_u = Malaria_prevalence + 1.96*sqrt(Malaria_prevalence * (1 - Malaria_prevalence) / Malaria_prev_n_tested))

df_species <- read_excel("data/data_RCT.xlsx", sheet = "site_params")

## read in the parameter efficacy estimates for nets
dat_res_pyr <- read.csv("parameters/pyrethroid_only_nets.csv",header=TRUE) 
dat_res_pbo <- read.csv("parameters/pyrethroid_pbo_nets.csv",header=TRUE) 
dat_res_pp <- read.csv("parameters/pyrethroid_pyrrole_nets.csv",header=TRUE)
dat_res_ppf <- read.csv("parameters/ppf_itn.csv") %>% mutate(resistance = 1 - bioassay_mortality)

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
# include a warm up period of ~30 years
sim_length <- 50 * 365

start_date <- as.Date("01/01/1980", format = "%d/%m/%Y")

# Benin times
# March 20th 2020 is the assumed intervention time - when nets are rolled out
# October to November 2019 - baseline time
mass_dist_dates_T <- c(as.Date("01/01/2015", format = "%d/%m/%Y"))
mass_dist_times_T <- as.vector(difftime(mass_dist_dates_T, start_date))

mass_dist_dates_B <- c(as.Date("01/01/2017", format = "%d/%m/%Y"))
mass_dist_times_B <- as.vector(difftime(mass_dist_dates_B, start_date))
  
intervention_date_B <- as.Date("20/03/2020", format = "%d/%m/%Y") # between March 19 and 22, 2020,LLINs were distributed
baseline_start_date_B <- as.Date("01/11/2019", format = "%d/%m/%Y") # baseline cross sectional survey was in October and November 2019 - calculated the mean prevalence on the baseline +/- 15 days

int_time_B <- difftime(intervention_date_B, start_date)[[1]]
baseline_time_B <- difftime(baseline_start_date_B, start_date)[[1]]

# Tanzania times

intervention_date_T <- as.Date("27/01/2019", format = "%d/%m/%Y") # LLINs were distributed among households between Jan 26 and Jan 28, 2019.
baseline_start_date_T <- as.Date("15/10/2018", format = "%d/%m/%Y") # we did cross-sectional surveys to collect data on household factors and malaria infection prevalence at baseline (October 2018) 

int_time_T <- difftime(intervention_date_T, start_date)[[1]]
baseline_time_T <- difftime(baseline_start_date_T, start_date)[[1]]

df$Time_months <- as.numeric(df$Time_months)

df <- rbind(subset(df, Country == "Tanzania") %>% rowwise() %>% mutate(date = if(is.na(Time_months) == 1){baseline_start_date_T}else{intervention_date_T + months(Time_months)}),
            subset(df, Country == "Benin") %>% rowwise() %>% mutate(date = if(is.na(Time_months) == 1){baseline_start_date_B}else{intervention_date_B + months(Time_months)})
            )
                                                  
df_B <- subset(df, Country == "Benin")

df_T <- subset(df, Country == "Tanzania")

# changing the nets each day
tu_diff_time_in <- 30

top_up_IG2_B <- get_top_up(net = "IG2",
                           df = df,
                           country = "Benin",
                           tu_diff_time = tu_diff_time_in,
                           sim_length = sim_length,
                           int_bed_net_time = int_time_B,
                           int_date = intervention_date_B)

top_up_RG_B <- get_top_up(net = "RG",
                          df = df,
                          country = "Benin",
                          tu_diff_time = tu_diff_time_in,
                          sim_length = sim_length,
                          int_bed_net_time = int_time_B,
                          int_date = intervention_date_B)

top_up_p_only_B <- get_top_up(net = "Pyrethroid_only",
                            df = df,
                            country = "Benin",
                            tu_diff_time = tu_diff_time_in,
                            sim_length = sim_length,
                            int_bed_net_time = int_time_B,
                            int_date = intervention_date_B)

top_up_IG2_T <- get_top_up(net = "IG2",
                           df = df,
                           country = "Tanzania",
                           tu_diff_time = tu_diff_time_in,
                           sim_length = sim_length,
                           int_bed_net_time = int_time_T,
                           int_date = intervention_date_T)

top_up_RG_T <- get_top_up(net = "RG",
                           df = df,
                           country = "Tanzania",
                           tu_diff_time = tu_diff_time_in,
                           sim_length = sim_length,
                           int_bed_net_time = int_time_T,
                          int_date = intervention_date_T)

top_up_PBO_T <- get_top_up(net = "PBO",
                           df = df,
                           country = "Tanzania",
                           tu_diff_time = tu_diff_time_in,
                           sim_length = sim_length,
                           int_bed_net_time = int_time_T,
                           int_date = intervention_date_T)

top_up_p_only_T <- get_top_up(net = "Pyrethroid_only",
                           df = df,
                           country = "Tanzania",
                           tu_diff_time = tu_diff_time_in,
                           sim_length = sim_length,
                           int_bed_net_time = int_time_T,
                           int_date = intervention_date_T)

top_up_all <- list("top_up_IG2_B" = top_up_IG2_B,
                 "top_up_RG_B" = top_up_RG_B,
                 "top_up_p_only_B" = top_up_p_only_B,
                 "top_up_IG2_T" = top_up_IG2_T,
                 "top_up_RG_T" = top_up_RG_T,
                 "top_up_PBO_T" = top_up_PBO_T,
                 "top_up_p_only_T" = top_up_p_only_T)

#saveRDS(top_up_all, file = "parameters/top_up_fits.rds")
top_up_all <- readRDS(file = "parameters/top_up_fits.rds")

top_up_IG2_B <- top_up_all$top_up_IG2_B
top_up_RG_B <- top_up_all$top_up_RG_B
top_up_p_only_B <- top_up_all$top_up_p_only_B
top_up_IG2_T <- top_up_all$top_up_IG2_T
top_up_RG_T <- top_up_all$top_up_RG_T
top_up_PBO_T <- top_up_all$top_up_PBO_T
top_up_p_only_T <- top_up_all$top_up_p_only_T

# plotting
png(file = "figures/net_cover.png", height = 700, width = 1100)
cowplot::plot_grid(
  top_up_p_only_T$plot + 
    ggtitle("Tanzania: pyrethroid-only") + 
    theme(legend.position = "none"),
  top_up_IG2_T$plot +
    ggtitle("Tanzania: pyrethroid-pyrrole") + 
    theme(legend.position = "none"),
  top_up_PBO_T$plot +
    ggtitle("Tanzania: pyrethroid-PBO") + 
    theme(legend.position = "none"),
  top_up_IG2_B$plot +
    ggtitle("Benin: pyrethroid-pyrrole") + 
    theme(legend.position = "none"),
  top_up_p_only_B$plot + 
    ggtitle("Benin: pyrethroid-only") +
    theme(legend.position = "none"),
  get_legend(top_up_IG2_T$plot),
  labels = c("A", "B", "C", "D", "E", "")
)
dev.off()

cowplot::plot_grid(
  cowplot::plot_grid(
  top_up_p_only_T$plot + 
    ggtitle("Tanzania: pyrethroid-only") + 
    theme(legend.position = "none"),
  top_up_RG_T$plot +
    ggtitle("Tanzania: pyrethroid-pyriproxyfen") + 
    theme(legend.position = "none"),
  top_up_p_only_B$plot +
    ggtitle("Benin: pyrethroid-only") + 
    theme(legend.position = "none"),
  top_up_RG_B$plot + 
    ggtitle("Benin: pyrethroid-pyriproxyfen") +
    theme(legend.position = "none")),
  get_legend(top_up_p_only_T$plot), ncol = 2, rel_widths = c(1, 0.2)
)

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
                             dat_res_pyr_in = dat_res_pyr,
                             
                             human_population = 10000,
                             year = 365,
                             sim_years = sim_length/365,
                             previous_mass_dist_times
                             ){
  
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
                                
                                # each row needs to show the efficacy parameter across years (and cols are diff mosquito)
                                # gambiae, coluzzi, arabiensis, funestus, coustani
                                # no resistance assumed for coustani
                                
                                dn0 = matrix(c(dn0, dn0, dn0,
                                               rep(dn0, n_mass_dist * 3)), nrow = (1 + n_mass_dist), ncol = 3),
                                rn = matrix(c(rn0, rn0, rn0,
                                              rep(rn0, n_mass_dist * 3)), nrow = (1 + n_mass_dist), ncol = 3),
                                rnm = matrix(c(rnm, rnm, rnm,
                                               rep(rnm, n_mass_dist * 3)), nrow = (1 + n_mass_dist), ncol = 3),
                                gamman = rep(gamman * 365,  (1 + n_mass_dist)) 
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

T_vals <- data.frame("Location_in" = rep("Misungwi", 4),
                     "Net_in" = c("IG2", "PBO", "RG", "Pyrethroid_only"),
                     "top_up_name" = c("top_up_IG2_T", "top_up_PBO_T", "top_up_RG_T", "top_up_p_only_T"),
                     "Country_in" = rep("Tanzania", 4))

cl <- makeCluster(4)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("T_site", "int_time_T", "baseline_time_T", "mass_dist_times_T",
                                   "df", "df_species_in", 
                                   "dat_res_pyr", "sim_length", 
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
    dat_res_pyr_in = dat_res_pyr,
    human_population = 10000,
    year = 365,
    sim_years = sim_length/365,
    previous_mass_dist_times = mass_dist_times_T)},
    error = function(cond){
      return(NA)
    })
}
saveRDS(start_EIR_T, 
        file = "data/start_EIR_T.rds")

stopCluster(cl)

start_EIR_T <- readRDS(file = "data/start_EIR_T.rds")

B_vals <- data.frame("Location_in" = c("Cove", "Zagnanado", "Ouinhi"),
                     "Net_in" = c("RG", "IG2", "Pyrethroid_only"),
                     "top_up_name" = c("top_up_RG_B", "top_up_IG2_B", "top_up_p_only_B"),
                     "Country_in" = rep("Benin", 3))

cl <- makeCluster(3)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("B_site", "int_time_B", "baseline_time_B", "mass_dist_times_B",
                                   "df", "df_species_in", 
                                   "dat_res_pyr", "sim_length", 
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
    intervention_time = int_time_B,
    baseline_start_time = baseline_time_B,
    df_in = as.data.frame(df),
    prop_species_in = as.data.frame(df_species_in),
    dat_res_pyr_in = dat_res_pyr,
    human_population = 10000,
    year = 365,
    sim_years = sim_length/365,
    previous_mass_dist_times = mass_dist_times_B)},
    error = function(cond){
      return(NA)
    })
}
saveRDS(start_EIR_B, 
        file = "data/start_EIR_B.rds")

stopCluster(cl)

start_EIR_B <- readRDS(file = "data/start_EIR_B.rds")

sim_forward <- function(start_EIR,
                        
                        Location_in,
                        Net_in,
                        Country_in,
                        
                        int_Net_in,
                        
                        top_up_name,
                        
                        top_up_net,
                        
                        sites_in,
                        
                        intervention_time,
                        baseline_pred_days,
                        
                        df_in = as.data.frame(df),
                        prop_species_in = as.data.frame(df_species_in),
                        
                        dat_res_pyr_in = dat_res_pyr,
                        dat_res_pp_in = dat_res_pp,
                        dat_res_pbo_in = dat_res_pbo,
                        dat_res_ppf_in = dat_res_ppf,
                        
                        human_population = 10000,
                        year = 365,
                        sim_years = sim_length/365,
                        em_diff = NULL,
                        bioassay_uncertainty = "middle",
                        
                        previous_mass_dist_times,
                        trial_net_cov_in = NA){
  
  print(paste0(Location_in, ": ", int_Net_in))
  
  top_up <- get(top_up_name)
  
  # specifying the retention parameter
  # if no top up then let the nets decline
  retention <- if(top_up_net == "none"){
    case_when(bioassay_uncertainty == "middle" ~ median(rstan::extract(top_up$fit_base, "retention")$retention) * 365,
              bioassay_uncertainty == "lower" ~ quantile(rstan::extract(top_up$fit_base, "retention")$retention, prob = 0.975)[[1]] * 365,
              bioassay_uncertainty == "upper" ~ quantile(rstan::extract(top_up$fit_base, "retention")$retention, prob = 0.025)[[1]] * 365)
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
  dn0_name <- case_when(bioassay_uncertainty == "middle" ~ "dn0_med",
                        bioassay_uncertainty == "lower" ~ "dn0_lo10",
                        bioassay_uncertainty == "upper" ~ "dn0_up90")
  
  # if bioassay mortality is lower then the lower repellency probability estimate is used else the higher value if used
  rn0_name <- case_when(bioassay_uncertainty == "middle" ~ "rn0_med",
                        bioassay_uncertainty == "lower" ~ "rn0_lo10",
                        bioassay_uncertainty == "upper" ~ "rn0_up90")
    
  # if bioassay mortality is lower then the longer bed net half life is used else the higher value if used
  gamman_name <- case_when(bioassay_uncertainty == "middle" ~ "gamman_med",
                           bioassay_uncertainty == "lower" ~ "gamman_up90",
                           bioassay_uncertainty == "upper" ~ "gamman_lo10")
  
  bioassay_index <- which(round(dat_res_pyr_in$bioassay_mortality, digits = 2) == bioassay_mortality)
  dn0 <- dat_res_pyr_in[bioassay_index, dn0_name]
  rn0 <- dat_res_pyr_in[bioassay_index, rn0_name]
  rnm <- 0.24
  gamman <- dat_res_pyr_in[bioassay_index, gamman_name] * 365
  
  dat_int <- if(int_Net_in == "IG2"){dat_res_pp_in} else if(int_Net_in == "PBO"){dat_res_pbo_in} else if(int_Net_in == "PPF"){dat_res_ppf_in} else if(int_Net_in == "Pyrethroid_only"){dat_res_pyr_in}
    
  int_bioassay_index <- which(round(dat_int$bioassay_mortality, digits = 2) == bioassay_mortality)
  dn0_int <- dat_int[int_bioassay_index, dn0_name]
  rn0_int <- dat_int[int_bioassay_index, rn0_name]
  gamman_int <- dat_int[int_bioassay_index, gamman_name] * 365
  
  dat_top_up <- if(top_up_net == "IG2"){dat_res_pp_in}else if(top_up_net == "PBO"){dat_res_pbo_in} else if(top_up_net == "PPF"){dat_res_ppf_in} else if(top_up_net == "Pyrethroid_only"){dat_res_pyr_in}
  
  if(top_up_net != "none"){
    tu_bioassay_index <- which(round(dat_top_up$bioassay_mortality, digits = 2) == bioassay_mortality)
    dn0_tu <- dat_top_up[tu_bioassay_index, dn0_name]
    rn0_tu <- dat_top_up[tu_bioassay_index, rn0_name]
    gamman_tu <- dat_top_up[tu_bioassay_index, gamman_name] * 365
  }
  
  
  bednetparams <- simparams
  
  init_cov <- min(df_in[baseline_index, "Bed_net_use_both"], top_up$mean_bed_net_use_both) # 
  covs_df <- top_up$covs_df
  
  n_mass_dist <- length(previous_mass_dist_times)
  
  if(top_up_net == "none"){
    
    trial_net_cov <- if(is.na(trial_net_cov_in) == FALSE){trial_net_cov_in}else{top_up$mean_bed_net_use_both}
    
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
                                  gamman = c(rep(gamman, n_init_times), gamman_int)
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
                      
  # less top up if bioassay mortality is higher
  top_up_name <- case_when(bioassay_uncertainty == "middle" ~ "top_up",
                           bioassay_uncertainty == "lower" ~ "top_up_l",
                           bioassay_uncertainty == "upper" ~ "top_up_u"
                           )
  
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

                                # dn0 = matrix(c(dn0, dn0_int, rep(dn0_tu, n_top_up), 
                                #                dn0, dn0_int, rep(dn0_tu, n_top_up),
                                #                dn0, dn0_int, rep(dn0_tu, n_top_up)), nrow = n_top_up + 2, ncol = 3),
                                # rn = matrix(c(rn0, rn0_int, rep(rn0_tu, n_top_up), 
                                #               rn0, rn0_int, rep(rn0_tu, n_top_up),
                                #               rn0, rn0_int, rep(rn0_tu, n_top_up)), nrow = n_top_up + 2, ncol = 3),
                                # rnm = matrix(c(rnm, rnm, rep(rnm, n_top_up), 
                                #                rnm, rnm, rep(rnm, n_top_up),
                                #                rnm, rnm, rep(rnm, n_top_up)), nrow = n_top_up + 2, ncol = 3),
                                
                                dn0 = matrix(c(dn0, dn0_mass_dist, dn0_vals,
                                               dn0, dn0_mass_dist, dn0_vals,
                                               dn0, dn0_mass_dist, dn0_vals), nrow = length(int_times) + n_mass_dist + 1, ncol = 3),
                                
                                rn = matrix(c(rn0, rn0_mass_dist, rn0_vals,
                                              rn0, rn0_mass_dist, rn0_vals,
                                              rn0, rn0_mass_dist, rn0_vals), nrow = length(int_times) + n_mass_dist + 1, ncol = 3),
                                
                                rnm = matrix(c(rnm, rnm_mass_dist, rep(rnm, length(int_times)),
                                               rnm, rnm_mass_dist, rep(rnm, length(int_times)),
                                               rnm, rnm_mass_dist, rep(rnm, length(int_times))), nrow = length(int_times) + n_mass_dist + 1, ncol = 3),

                                gamman = c(gamman, rep(gamman, n_mass_dist), gamman_in))
    
  }
  
  if(is.null(em_diff) != 1){
    bednetparams_1 <- set_habitat_management(
      parameters = bednetparams_1,
      habitat_management_timesteps = c(intervention_time), # TRIAL nets first introduced on year 4
      lsm_new_eqm = matrix(rep(em_diff, 3), nrow = 1, ncol = 3), # assumes
      lsm_rate_alpha = matrix(rep(-4, 3), nrow = 1, ncol = 3),
      lsm_rate_beta = matrix(rep(0.1, 3), nrow = 1, ncol = 3)
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
  
  simparams_eq <- set_equilibrium(parameters = bednetparams_1, init_EIR = start_EIR)
  out <- run_simulation(sim_length, simparams_eq, correlations = correlationsb1)
  return(out)
}

T_vals$start_EIR <- unlist(start_EIR_T)

# running the simulations
T_vals_pred_u <- rbind(T_vals %>% mutate(int_Net_in = "IG2", start_EIR = unlist(start_EIR_T)),
                     T_vals %>% mutate(int_Net_in = "PBO", start_EIR = unlist(start_EIR_T)),
                     T_vals %>% mutate(int_Net_in = "Pyrethroid_only", start_EIR = unlist(start_EIR_T)))

T_vals_pred_u <- subset(T_vals_pred_u, Net_in != "RG")

T_vals_pred_b <- rbind(T_vals_pred_u %>% mutate(bioassay_uncertainty = "middle"),
                       T_vals_pred_u %>% mutate(bioassay_uncertainty = "lower"),
                       T_vals_pred_u %>% mutate(bioassay_uncertainty = "upper")) %>% as.data.frame()

T_vals_pred <- bind_rows(
  lapply(seq(1, 3, 1), function(i, T_vals_pred_b){
    T_vals_pred_b %>% mutate(rep = i)}, 
    T_vals_pred_b = T_vals_pred_b)
  )

cl <- makeCluster(5)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("T_site", "int_time_T", "baseline_time_T", "df", "df_species_in", 
                                   "dat_res_pyr", "sim_length", "dat_res_pp", "dat_res_pbo",
                                   "top_up_IG2_T", "top_up_PBO_T", "mass_dist_times_T",
                                   "top_up_RG_T", "top_up_p_only_T",
                                   "decline_d0", "decline_r0",
                                   "calc_base_mean"))
pred_prev_T <- foreach(i=1:nrow(T_vals_pred),
                       .packages = (.packages())
) %dopar% {
  tryCatch({sim_forward(start_EIR = T_vals_pred[i, "start_EIR"],
                        Location_in = T_vals_pred[i, "Location_in"],
                        Net_in = T_vals_pred[i, "Net_in"],
                        Country_in = T_vals_pred[i, "Country_in"],
                        int_Net_in = T_vals_pred[i, "int_Net_in"],
                        top_up_name = T_vals_pred[i, "top_up_name"],
                        top_up_net = "Pyrethroid_only",
                        sites_in = T_site,
                        intervention_time = int_time_T,
                        baseline_pred_days = baseline_time_T,
                        bioassay_uncertainty = T_vals_pred[i, "bioassay_uncertainty"],
                        previous_mass_dist_times = mass_dist_times_T)
  },
  error = function(cond){
    return(NA)
  })
}
saveRDS(pred_prev_T, 
        file = "data/pred_prev_T.rds")

stopCluster(cl)

pred_prev_T <- readRDS(file = "data/pred_prev_T.rds")

# top up simulations

# make differences in the parameter values - killing really extreme and then check that there are differences.

T_vals_tu <- T_vals %>% mutate(int_Net_in = Net_in, 
                  start_EIR = unlist(start_EIR_T),
                  top_up_net = Net_in)

T_vals_tu_s <- subset(T_vals_tu, Net_in != "RG")

T_vals_pred_tu_s <- rbind(T_vals_tu_s %>% mutate(bioassay_uncertainty = "middle"),
                       T_vals_tu_s %>% mutate(bioassay_uncertainty = "lower"),
                       T_vals_tu_s %>% mutate(bioassay_uncertainty = "upper")) %>% as.data.frame()

T_vals_pred_tu <- bind_rows(
  lapply(seq(1, 3, 1), function(i, T_vals_pred_b){
    T_vals_pred_b %>% mutate(rep = i)}, 
    T_vals_pred_b = T_vals_pred_tu_s)
)

cl <- makeCluster(5)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("T_site", "int_time_T", "baseline_time_T",  "mass_dist_times_T",
                                   "df", "df_species_in", 
                                   "dat_res_pyr", "sim_length", "dat_res_pp", "dat_res_pbo",
                                   "top_up_IG2_T", "top_up_PBO_T", 
                                   "top_up_RG_T", "top_up_p_only_T",
                                   "decline_d0", "decline_r0",
                                   "calc_base_mean"))
pred_prev_T_tu <- foreach(i=1:nrow(T_vals_pred_tu),
                       .packages = (.packages())
) %dopar% {
  tryCatch({sim_forward(start_EIR = T_vals_pred_tu[i, "start_EIR"],
                        Location_in = T_vals_pred_tu[i, "Location_in"],
                        Net_in = T_vals_pred_tu[i, "Net_in"],
                        Country_in = T_vals_pred_tu[i, "Country_in"],
                        int_Net_in = T_vals_pred_tu[i, "int_Net_in"],
                        top_up_name = T_vals_pred_tu[i, "top_up_name"],
                        top_up_net = T_vals_pred_tu[i, "top_up_net"],
                        sites_in = T_site,
                        intervention_time = int_time_T,
                        baseline_pred_days = baseline_time_T,
                        bioassay_uncertainty = T_vals_pred_tu[i, "bioassay_uncertainty"],
                        previous_mass_dist_times = mass_dist_times_T)
  },
  error = function(cond){
    return(NA)
  })
}
saveRDS(pred_prev_T_tu, 
        file = "data/pred_prev_T_tu.rds")

stopCluster(cl)

pred_prev_T_tu <- readRDS(file = "data/pred_prev_T_tu.rds")

# model with no top up
T_vals_pred_tu_none <- T_vals_pred_tu %>% mutate(top_up_net = "none",
                                                 trial_net_cov_in = 0)

cl <- makeCluster(5)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("T_site", "int_time_T", "baseline_time_T",  "mass_dist_times_T",
                                   "df", "df_species_in", 
                                   "dat_res_pyr", "sim_length", "dat_res_pp", "dat_res_pbo",
                                   "top_up_IG2_T", "top_up_PBO_T", 
                                   "top_up_RG_T", "top_up_p_only_T",
                                   "decline_d0", "decline_r0",
                                   "calc_base_mean"))
pred_prev_T_tu_none <- foreach(i=1:nrow(T_vals_pred_tu_none),
                          .packages = (.packages())
) %dopar% {
  tryCatch({sim_forward(start_EIR = T_vals_pred_tu_none[i, "start_EIR"],
                        Location_in = T_vals_pred_tu_none[i, "Location_in"],
                        Net_in = T_vals_pred_tu_none[i, "Net_in"],
                        Country_in = T_vals_pred_tu_none[i, "Country_in"],
                        int_Net_in = T_vals_pred_tu_none[i, "int_Net_in"],
                        top_up_name = T_vals_pred_tu_none[i, "top_up_name"],
                        top_up_net = T_vals_pred_tu_none[i, "top_up_net"],
                        sites_in = T_site,
                        intervention_time = int_time_T,
                        baseline_pred_days = baseline_time_T,
                        bioassay_uncertainty = T_vals_pred_tu_none[i, "bioassay_uncertainty"],
                        previous_mass_dist_times = mass_dist_times_T,
                        trial_net_cov_in = T_vals_pred_tu_none[i, "trial_net_cov_in"])
  },
  error = function(cond){
    return(NA)
  })
}
saveRDS(pred_prev_T_tu_none, 
        file = "data/pred_prev_T_tu_none.rds")

stopCluster(cl)

pred_prev_T_tu_none <- readRDS(file = "data/pred_prev_T_tu_none.rds")

# Benin
B_vals_pred_u <- rbind(B_vals %>% mutate(int_Net_in = "IG2", start_EIR = unlist(start_EIR_B)),
                     B_vals %>% mutate(int_Net_in = "PBO", start_EIR = unlist(start_EIR_B)),
                     B_vals %>% mutate(int_Net_in = "Pyrethroid_only", start_EIR = unlist(start_EIR_B)))

B_vals_pred_u <- subset(B_vals_pred_u, Net_in != "RG")

B_vals_pred_b <- rbind(B_vals_pred_u %>% mutate(bioassay_uncertainty = "middle"),
                       B_vals_pred_u %>% mutate(bioassay_uncertainty = "lower"),
                       B_vals_pred_u %>% mutate(bioassay_uncertainty = "upper")) 

B_vals_pred <- bind_rows(
  lapply(seq(1, 3, 1), function(i, B_vals_pred_b){
    B_vals_pred_b %>% mutate(rep = i)}, 
    B_vals_pred_b = B_vals_pred_b)
)

cl <- makeCluster(5)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("B_site", "int_time_B", "baseline_time_B",  "mass_dist_times_B",
                                   "df", "df_species_in", 
                                   "dat_res_pyr", "sim_length", "dat_res_pp", "dat_res_pbo",
                                   "top_up_IG2_B",
                                   "top_up_RG_B", "top_up_p_only_B",
                                   "decline_d0", "decline_r0",
                                   "calc_base_mean"))
pred_prev_B <- foreach(i=1:nrow(B_vals_pred),
                       .packages = (.packages())
) %dopar% {
  tryCatch({sim_forward(start_EIR = B_vals_pred[i, "start_EIR"],
                        Location_in = B_vals_pred[i, "Location_in"],
                        Net_in = B_vals_pred[i, "Net_in"],
                        Country_in = B_vals_pred[i, "Country_in"],
                        int_Net_in = B_vals_pred[i, "int_Net_in"],
                        top_up_name = B_vals_pred[i, "top_up_name"],
                        top_up_net = "Pyrethroid_only",
                        sites_in = B_site,
                        intervention_time = int_time_B,
                        baseline_pred_days = baseline_time_B,
                        bioassay_uncertainty = B_vals_pred[i, "bioassay_uncertainty"],
                        previous_mass_dist_times = mass_dist_times_B)
  },
  error = function(cond){
    return(NA)
  })
}
saveRDS(pred_prev_B, 
        file = "data/pred_prev_cf_B.rds")

stopCluster(cl)

pred_prev_B <- readRDS(file = "data/pred_prev_cf_B.rds")

# top up simulations
B_vals_tu <- B_vals %>% mutate(int_Net_in = Net_in, 
                               start_EIR = unlist(start_EIR_B),
                               top_up_net = Net_in)

B_vals_tu_s <- subset(B_vals_tu, Net_in != "RG")

B_vals_pred_tu_s <- rbind(B_vals_tu_s %>% mutate(bioassay_uncertainty = "middle"),
                          B_vals_tu_s %>% mutate(bioassay_uncertainty = "lower"),
                          B_vals_tu_s %>% mutate(bioassay_uncertainty = "upper")) %>% as.data.frame()

B_vals_pred_tu <- bind_rows(
  lapply(seq(1, 3, 1), function(i, T_vals_pred_b){
    T_vals_pred_b %>% mutate(rep = i)}, 
    T_vals_pred_b = B_vals_pred_tu_s)
)

cl <- makeCluster(5)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("B_site", "int_time_B", "baseline_time_B",  "mass_dist_times_B",
                                   "df", "df_species_in", 
                                   "dat_res_pyr", "sim_length", "dat_res_pp", "dat_res_pbo",
                                   "top_up_IG2_B",
                                   "top_up_RG_B", "top_up_p_only_B",
                                   "decline_d0", "decline_r0",
                                   "calc_base_mean"))
pred_prev_B_tu <- foreach(i=1:nrow(B_vals_pred_tu),
                       .packages = (.packages())
) %dopar% {
  tryCatch({sim_forward(start_EIR = B_vals_pred_tu[i, "start_EIR"],
                        Location_in = B_vals_pred_tu[i, "Location_in"],
                        Net_in = B_vals_pred_tu[i, "Net_in"],
                        Country_in = B_vals_pred_tu[i, "Country_in"],
                        int_Net_in = B_vals_pred_tu[i, "int_Net_in"],
                        top_up_name = B_vals_pred_tu[i, "top_up_name"],
                        top_up_net = B_vals_pred_tu[i, "top_up_net"],
                        sites_in = B_site,
                        intervention_time = int_time_B,
                        baseline_pred_days = baseline_time_B,
                        bioassay_uncertainty = B_vals_pred_tu[i, "bioassay_uncertainty"],
                        previous_mass_dist_times = mass_dist_times_B)
  },
  error = function(cond){
    return(NA)
  })
}
saveRDS(pred_prev_B_tu, 
        file = "data/pred_prev_B_tu.rds")

stopCluster(cl)

pred_prev_B_tu <- readRDS(file = "data/pred_prev_B_tu.rds")

B_vals_pred_tu_none <- B_vals_pred_tu %>% mutate(top_up_net = "none",
                                                 trial_net_cov_in = 0)
                             
cl <- makeCluster(5)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("B_site", "int_time_B", "baseline_time_B", "mass_dist_times_B",
                                   "df", "df_species_in", 
                                   "dat_res_pyr", "sim_length", "dat_res_pp", "dat_res_pbo",
                                   "top_up_IG2_B",
                                   "top_up_RG_B", "top_up_p_only_B",
                                   "decline_d0", "decline_r0",
                                   "calc_base_mean"))
pred_prev_B_tu_none <- foreach(i=1:nrow(B_vals_pred_tu_none),
                          .packages = (.packages())
) %dopar% {
  tryCatch({sim_forward(start_EIR = B_vals_pred_tu_none[i, "start_EIR"],
                        Location_in = B_vals_pred_tu_none[i, "Location_in"],
                        Net_in = B_vals_pred_tu_none[i, "Net_in"],
                        Country_in = B_vals_pred_tu_none[i, "Country_in"],
                        int_Net_in = B_vals_pred_tu_none[i, "int_Net_in"],
                        top_up_name = B_vals_pred_tu_none[i, "top_up_name"],
                        top_up_net = B_vals_pred_tu_none[i, "top_up_net"],
                        sites_in = B_site,
                        intervention_time = int_time_B,
                        baseline_pred_days = baseline_time_B,
                        bioassay_uncertainty = B_vals_pred_tu_none[i, "bioassay_uncertainty"],
                        previous_mass_dist_times = mass_dist_times_B,
                        trial_net_cov_in = B_vals_pred_tu_none[i, "trial_net_cov_in"])
  },
  error = function(cond){
    return(NA)
  })
}
saveRDS(pred_prev_B_tu_none, 
        file = "data/pred_prev_B_tu_none.rds")

stopCluster(cl)

pred_prev_B_tu_none <- readRDS(file = "data/pred_prev_B_tu_none.rds")

# RG model
devtools::dev_mode(on = TRUE, path = getOption("devtools.path"))
detach("package:malariasimulation", unload=TRUE)
library(malariasimulation); library(malariaEquilibrium)

em_mos <- readxl::read_excel("data/adult_reduction_benin.xlsx", sheet = 1)
mean_M_decline_all <- mean(em_mos$PPF_no_usage_decline / em_mos$PYR_no_usage_decline)

mean_M_decline_arm <- mean((em_mos$PYR_no_usage_decline - (em_mos$PYR_no_usage_decline - em_mos$PPF_no_usage_decline)*0.57)/em_mos$PYR_no_usage_decline)

out_ppf <- sim_forward(start_EIR = start_EIR_B[[1]],
                       Location_in = "Cove",
                       Net_in = "RG",
                       Country_in = "Benin",
                       sites_in = B_site,
                       int_Net_in = "RG",
                       top_up_name = "top_up_RG_B",
                       top_up_net = "pyrethroid_only",
                       intervention_time = int_time_B,
                       baseline_pred_days = baseline_time_B,
                       em_diff = mean_M_decline_arm)

B_vals_RG <- B_vals %>% mutate(int_Net_in = "RG", start_EIR = unlist(start_EIR_B))

pred_prev_B_RG <- foreach(i=1:nrow(B_vals_RG),
                       .packages = (.packages())
) %do% {sim_forward(start_EIR = B_vals_RG[i, "start_EIR"],
                        Location_in = B_vals_RG[i, "Location_in"],
                        Net_in = B_vals_RG[i, "Net_in"],
                        Country_in = B_vals_RG[i, "Country_in"],
                        int_Net_in = B_vals_RG[i, "int_Net_in"],
                        top_up_name = B_vals_RG[i, "top_up_name"],
                        top_up_net = "Pyrethroid_only",
                        sites_in = B_site,
                        intervention_time = int_time_B,
                        baseline_pred_days = baseline_time_B,
                        em_diff = mean_M_decline_arm)
}

saveRDS(pred_prev_B_RG, 
        file = "data/pred_prev_cf_RG_B.rds")

devtools::dev_mode(on = FALSE, path = getOption("devtools.path"))

pred_prev_B_RG <- readRDS(file = "data/pred_prev_cf_RG_B.rds")

# formatting the results
# functions to extract the results
summary_pfpr_0_100_all <- function(output){
  return(output[,'n_detect_0_36500'] / output[,'n_0_36500'])
}

summary_pfpr_0.5_5 <- function(output){
  return(output[,'n_detect_182.5_5110'] / output[,'n_182.5_5110'])
}

for(i in 1:length(pred_prev_B)){
  pred_prev_B[[i]]$date <- seq(start_date, by = "day", length.out = length(pred_prev_B[[i]]$timestep))
}

for(i in 1:length(pred_prev_B_tu)){
  pred_prev_B_tu[[i]]$date <- seq(start_date, by = "day", length.out = length(pred_prev_B_tu[[i]]$timestep))
}

for(i in 1:length(pred_prev_B_tu_none)){
  pred_prev_B_tu_none[[i]]$date <- seq(start_date, by = "day", length.out = length(pred_prev_B_tu_none[[i]]$timestep))
}

for(i in 1:length(pred_prev_T)){
  pred_prev_T[[i]]$date <- seq(start_date, by = "day", length.out = length(pred_prev_T[[1]]$timestep))
}

for(i in 1:length(pred_prev_T_tu)){
  pred_prev_T_tu[[i]]$date <- seq(start_date, by = "day", length.out = length(pred_prev_T_tu[[1]]$timestep))
}

for(i in 1:length(pred_prev_T_tu_none)){
  pred_prev_T_tu_none[[i]]$date <- seq(start_date, by = "day", length.out = length(pred_prev_T_tu_none[[1]]$timestep))
}

extract_prev <- function(Net_in, 
                              int_Net_in,
                              Country_in,
                              bioassay_uncertainty,
                              vals_pred, 
                              pred_prev,
                              top_up_net = NA,
                              n_rep = 3,
                         trial_net_cov_in = NULL){
  
  #p <- if(bioassay_uncertainty == "middle"){0.5} else if(bioassay_uncertainty == "lower"){0.025} else{0.975}
  
  if(is.null(trial_net_cov_in) == 1){
  inds <- if(is.na(top_up_net) == 1){
    which(vals_pred$Net_in == Net_in & 
          vals_pred$int_Net_in == int_Net_in &
          vals_pred$bioassay_uncertainty == bioassay_uncertainty)
  } else{
    which(vals_pred$Net_in == Net_in & 
            vals_pred$int_Net_in == int_Net_in &
            vals_pred$bioassay_uncertainty == bioassay_uncertainty &
            vals_pred$top_up_net == top_up_net)
  }} else{
    
    inds <- if(is.na(top_up_net) == 1){
      which(vals_pred$Net_in == Net_in & 
              vals_pred$int_Net_in == int_Net_in &
              vals_pred$bioassay_uncertainty == bioassay_uncertainty &
              is.na(vals_pred$trial_net_cov_in) == is.na(trial_net_cov_in))
    } else{
      which(vals_pred$Net_in == Net_in & 
              vals_pred$int_Net_in == int_Net_in &
              vals_pred$bioassay_uncertainty == bioassay_uncertainty &
              vals_pred$top_up_net == top_up_net &
      is.na(vals_pred$trial_net_cov_in) == is.na(trial_net_cov_in))
    }
    
    
  }
  
  
  if(length(inds) != n_rep){return(NA)}else{
    
    prev_mat <- if(Country_in == "Tanzania"){sapply(inds, function(i, p){summary_pfpr_0.5_5(p[[i]])}, p = pred_prev)
    } else if(Country_in == "Benin"){sapply(inds, function(i, p){summary_pfpr_0_100_all(p[[i]])}, p = pred_prev)
    }
    
    return(apply(prev_mat, 1, mean))
  }
}

extract_mean_prev <- function(Net_in, 
                              int_Net_in,
                              Country_in,
                              vals_pred, 
                              pred_prev,
                              top_up_net = NA,
                              n_rep = 3,
                              trial_net_cov_in = NULL){
  
  
  inds <- if(is.null(trial_net_cov_in) == 1){if(is.na(top_up_net) == 1){
    which(vals_pred$Net_in == Net_in & 
            vals_pred$int_Net_in == int_Net_in)
  } else{
    which(vals_pred$Net_in == Net_in & 
            vals_pred$int_Net_in == int_Net_in &
            vals_pred$top_up_net == top_up_net)
  }} else{
    if(is.na(top_up_net) == 1){
      which(vals_pred$Net_in == Net_in & 
              vals_pred$int_Net_in == int_Net_in &
              is.na(vals_pred$trial_net_cov_in) == is.na(trial_net_cov_in))
    } else{
      which(vals_pred$Net_in == Net_in & 
              vals_pred$int_Net_in == int_Net_in &
              vals_pred$top_up_net == top_up_net &
              is.na(vals_pred$trial_net_cov_in) == is.na(trial_net_cov_in))
    }
  }
  
  if(length(inds) != n_rep*3){return(NA)}else{
    
    prev_mat <- if(Country_in == "Tanzania"){sapply(inds, function(i, p){summary_pfpr_0.5_5(p[[i]])}, p = pred_prev)
    } else if(Country_in == "Benin"){sapply(inds, function(i, p){summary_pfpr_0_100_all(p[[i]])}, p = pred_prev)
    }
    
    out <- as.data.frame(prev_mat) %>% mutate(mean = rowMeans(prev_mat),
                                              t = pred_prev[[1]]$timestep,
                                              date = pred_prev[[1]]$date) %>% 
      pivot_longer(!c(t, date), names_to = "rep", values_to = "prev") %>% 
      as.data.frame()
    
    return(out)
  }
}

arm_sims_T <- data.frame("t" = rep(pred_prev_T[[1]]$timestep, 3),
                            "date" = rep(pred_prev_T[[1]]$date, 3),
                            "net" = c(rep("IG2", sim_length),
                                      rep("PBO", sim_length),
                                      rep("Pyrethroid_only", sim_length)),
                            
                            "prev" = c(extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                                    Country_in = "Tanzania", bioassay_uncertainty = "middle",
                                                    vals_pred = T_vals_pred, pred_prev = pred_prev_T),
                                       
                                       extract_prev(Net_in = "PBO", int_Net_in = "PBO",
                                                    Country_in = "Tanzania", bioassay_uncertainty = "middle",
                                                    vals_pred = T_vals_pred, pred_prev = pred_prev_T),
                                       
                                       extract_prev(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only",
                                                    Country_in = "Tanzania", bioassay_uncertainty = "middle",
                                                    vals_pred = T_vals_pred, pred_prev = pred_prev_T)
                            ),
                            
                            "lower" = c(extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                                     Country_in = "Tanzania", bioassay_uncertainty = "lower",
                                                     vals_pred = T_vals_pred, pred_prev = pred_prev_T),
                                        
                                        extract_prev(Net_in = "PBO", int_Net_in = "PBO",
                                                     Country_in = "Tanzania", bioassay_uncertainty = "lower",
                                                     vals_pred = T_vals_pred, pred_prev = pred_prev_T),
                                        
                                        extract_prev(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only",
                                                     Country_in = "Tanzania", bioassay_uncertainty = "lower",
                                                     vals_pred = T_vals_pred, pred_prev = pred_prev_T)
                            ),
                            
                            "upper" = c(extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                                     Country_in = "Tanzania", bioassay_uncertainty = "upper",
                                                     vals_pred = T_vals_pred, pred_prev = pred_prev_T),
                                        
                                        extract_prev(Net_in = "PBO", int_Net_in = "PBO",
                                                     Country_in = "Tanzania", bioassay_uncertainty = "upper",
                                                     vals_pred = T_vals_pred, pred_prev = pred_prev_T),
                                        
                                        extract_prev(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only",
                                                     Country_in = "Tanzania", bioassay_uncertainty = "upper",
                                                     vals_pred = T_vals_pred, pred_prev = pred_prev_T)
                            )
)

arm_sims_B <- data.frame("t" = rep(pred_prev_B[[1]]$timestep, 2),
                            "date" = rep(pred_prev_B[[1]]$date, 2),
                            "net" = c(rep("IG2", sim_length),
                                      rep("Pyrethroid_only", sim_length)),
                            "prev" = c(extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                                    Country_in = "Benin", bioassay_uncertainty = "middle",
                                                    vals_pred = B_vals_pred, pred_prev = pred_prev_B),
                                       extract_prev(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only",
                                                    Country_in = "Benin", bioassay_uncertainty = "middle",
                                                    vals_pred = B_vals_pred, pred_prev = pred_prev_B)
                            ),
                            "lower" = c(extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                                     Country_in = "Benin", bioassay_uncertainty = "lower",
                                                     vals_pred = B_vals_pred, pred_prev = pred_prev_B),
                                        extract_prev(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only",
                                                     Country_in = "Benin", bioassay_uncertainty = "lower",
                                                     vals_pred = B_vals_pred, pred_prev = pred_prev_B)
                            ),
                            "upper" = c(extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                                     Country_in = "Benin", bioassay_uncertainty = "upper",
                                                     vals_pred = B_vals_pred, pred_prev = pred_prev_B),
                                        extract_prev(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only",
                                                     Country_in = "Benin", bioassay_uncertainty = "upper",
                                                     vals_pred = B_vals_pred, pred_prev = pred_prev_B)
                            )
)

# arm_sims_T <- rbind(extract_prev(Net_in = "IG2", int_Net_in = "IG2", 
#                                                       Country_in = "Tanzania",
#                                                       vals_pred = T_vals_pred, 
#                                                       pred_prev = pred_prev_T) %>% mutate(net = "IG2"),
#                                     extract_mean_prev(Net_in = "PBO", int_Net_in = "PBO", 
#                                                       Country_in = "Tanzania",
#                                                       vals_pred = T_vals_pred, 
#                                                       pred_prev = pred_prev_T) %>% mutate(net = "PBO"),
#                                     extract_mean_prev(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only", 
#                                                       Country_in = "Tanzania",
#                                                       vals_pred = T_vals_pred, 
#                                                       pred_prev = pred_prev_T) %>% mutate(net = "Pyrethroid_only")
#                                     )

# arm_sims_B <- rbind(extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
#                                       Country_in = "Benin",
#                                       vals_pred = B_vals_pred, 
#                                       pred_prev = pred_prev_B) %>% mutate(net = "IG2"),
#                     extract_mean_prev(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only", 
#                                       Country_in = "Benin",
#                                       vals_pred = B_vals_pred, 
#                                       pred_prev = pred_prev_B) %>% mutate(net = "Pyrethroid_only")
# )

# Tanzania
df_T_prev <- na.omit(df_T[,c("Net", "date", "Malaria_prevalence", "Malaria_prevalence_l", "Malaria_prevalence_u")])

T_arm_plot <- ggplot() +
  geom_vline(xintercept = intervention_date_T, linetype = 2, linewidth = 1, alpha = 0.5) +
  geom_ribbon(data = arm_sims_T,
              aes(x = date, y = prev, ymin = upper, ymax = lower, 
                 group = net, fill = net), alpha = 0.4) +
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

df_B_prev <- na.omit(df_B[,c("Net", "date", "Malaria_prevalence", "Malaria_prevalence_l", "Malaria_prevalence_u", "Location")])

B_arm_plot <- ggplot() +
  geom_vline(xintercept = intervention_date_B, linetype = 2, linewidth = 1, alpha = 0.5) +
  geom_ribbon(data = arm_sims_B, 
              aes(x = date, ymin = upper, ymax = lower, fill = net, group = net), 
              alpha = 0.4) +
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

plot_grid(T_arm_plot, B_arm_plot)

# pyrethroid-pyrrole arms with pyrethroid-pyrrole top up nets

# rbind(extract_prev(Net_in = "IG2", int_Net_in = "IG2", 
#                    Country_in = "Tanzania",
#                    vals_pred = T_vals_pred,
#                    
#                    pred_prev = pred_prev_T) %>% mutate(net = "IG2"),
#       
#       extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
#                         Country_in = "Tanzania",
#                         vals_pred = T_vals_pred_tu, 
#                         pred_prev = pred_prev_T_tu) %>% mutate(net = "IG2_tu"),
#       
#       extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
#                         Country_in = "Tanzania",
#                         vals_pred = T_vals_pred_tu_none, 
#                         pred_prev = pred_prev_T_tu_none,
#                         trial_net_cov_in = NA) %>% mutate(net = "IG2_tu_none"),
#       
#       extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
#                         Country_in = "Tanzania",
#                         vals_pred = T_vals_pred_tu_none, 
#                         pred_prev = pred_prev_T_tu_none,
#                         trial_net_cov_in = 0) %>% mutate(net = "None")
# )

T_tu_df <- data.frame("t" = rep(pred_prev_T[[1]]$timestep, 3),
                                    "date" = rep(pred_prev_T[[1]]$date, 3),
                                    "net" = c(rep("IG2", sim_length),
                                              rep("IG2_tu", sim_length),
                                              rep("None", sim_length)),
                                    
                                    "prev" = c(extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                                            Country_in = "Tanzania", bioassay_uncertainty = "middle",
                                                            vals_pred = T_vals_pred, pred_prev = pred_prev_T),
                                               
                                               extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                                            Country_in = "Tanzania", bioassay_uncertainty = "middle",
                                                            vals_pred = T_vals_pred_tu, pred_prev = pred_prev_T_tu),
                                               
                                               extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                                            Country_in = "Tanzania", bioassay_uncertainty = "middle",
                                                            vals_pred = T_vals_pred_tu_none, pred_prev = pred_prev_T_tu_none,
                                                            trial_net_cov_in = 0)
                                    ),
                                    
                                    "lower" = c(extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                                             Country_in = "Tanzania", bioassay_uncertainty = "lower",
                                                             vals_pred = T_vals_pred, pred_prev = pred_prev_T),
                                                
                                                extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                                             Country_in = "Tanzania", bioassay_uncertainty = "lower",
                                                             vals_pred = T_vals_pred_tu, pred_prev = pred_prev_T_tu),
                                                
                                                extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                                             Country_in = "Tanzania", bioassay_uncertainty = "lower",
                                                             vals_pred = T_vals_pred_tu_none, pred_prev = pred_prev_T_tu_none,
                                                             trial_net_cov_in = 0)
                                    ),
                                    
                                    "upper" = c(extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                                             Country_in = "Tanzania", bioassay_uncertainty = "upper",
                                                             vals_pred = T_vals_pred, pred_prev = pred_prev_T),
                                                
                                                extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                                             Country_in = "Tanzania", bioassay_uncertainty = "upper",
                                                             vals_pred = T_vals_pred_tu, pred_prev = pred_prev_T_tu),
                                                
                                                extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                                             Country_in = "Tanzania", bioassay_uncertainty = "upper",
                                                             vals_pred = T_vals_pred_tu_none, pred_prev = pred_prev_T_tu_none,
                                                             trial_net_cov_in = 0)
                                    )
)

B_tu_df <- data.frame("t" = rep(pred_prev_B[[1]]$timestep, 3),
                      "date" = rep(pred_prev_B[[1]]$date, 3),
                      "net" = c(rep("IG2", sim_length),
                                rep("IG2_tu", sim_length),
                                rep("None", sim_length)),
                      
                      "prev" = c(extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                              Country_in = "Benin", bioassay_uncertainty = "middle",
                                              vals_pred = B_vals_pred, pred_prev = pred_prev_B),
                                 
                                 extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                              Country_in = "Benin", bioassay_uncertainty = "middle",
                                              vals_pred = B_vals_pred_tu, pred_prev = pred_prev_B_tu),
                                 
                                 extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                              Country_in = "Benin", bioassay_uncertainty = "middle",
                                              vals_pred = B_vals_pred_tu_none, pred_prev = pred_prev_B_tu_none,
                                              trial_net_cov_in = 0)
                      ),
                      
                      "lower" = c(extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                               Country_in = "Benin", bioassay_uncertainty = "lower",
                                               vals_pred = B_vals_pred, pred_prev = pred_prev_B),
                                  
                                  extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                               Country_in = "Benin", bioassay_uncertainty = "lower",
                                               vals_pred = B_vals_pred_tu, pred_prev = pred_prev_B_tu),
                                  
                                  extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                               Country_in = "Benin", bioassay_uncertainty = "lower",
                                               vals_pred = B_vals_pred_tu_none, pred_prev = pred_prev_B_tu_none,
                                               trial_net_cov_in = 0)
                      ),
                      
                      "upper" = c(extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                               Country_in = "Benin", bioassay_uncertainty = "upper",
                                               vals_pred = B_vals_pred, pred_prev = pred_prev_B),
                                  
                                  extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                               Country_in = "Benin", bioassay_uncertainty = "upper",
                                               vals_pred = B_vals_pred_tu, pred_prev = pred_prev_B_tu),
                                  
                                  extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                               Country_in = "Benin", bioassay_uncertainty = "upper",
                                               vals_pred = B_vals_pred_tu_none, pred_prev = pred_prev_B_tu_none,
                                               trial_net_cov_in = 0)
                      )
)

# rbind(extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
#                                    Country_in = "Benin",
#                                    vals_pred = B_vals_pred, 
#                                    pred_prev = pred_prev_B) %>% mutate(net = "IG2"),
#                  
#                  extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
#                                    Country_in = "Benin",
#                                    vals_pred = B_vals_pred_tu, 
#                                    pred_prev = pred_prev_B_tu) %>% mutate(net = "IG2_tu"),
#                  
#                  extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
#                                    Country_in = "Benin",
#                                    vals_pred = B_vals_pred_tu_none, 
#                                    pred_prev = pred_prev_B_tu_none,
#                                    trial_net_cov_in = NA) %>% mutate(net = "IG2_tu_none"),
#                  
#                  extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
#                                    Country_in = "Benin",
#                                    vals_pred = B_vals_pred_tu_none, 
#                                    pred_prev = pred_prev_B_tu_none,
#                                    trial_net_cov_in = 0) %>% mutate(net = "None")
# )

cf_plot_T <- ggplot() +
  geom_ribbon(data = subset(T_tu_df, date <= intervention_date_T),
            aes(x = date, ymin = upper, ymax = lower, fill = net, group = net),
            alpha = 0.1) +
  
  geom_line(data = subset(T_tu_df, date <= intervention_date_T),
            aes(x = date, y = prev, col = net, group = net),
            alpha = 0.1, linewidth = 0.2) +
  
  geom_ribbon(data = subset(T_tu_df, date >= intervention_date_T),
            aes(x = date, ymin = upper, ymax = lower, fill = net, group = net),
            alpha = 0.4) +
  
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
                                 "Trial nets replaced\nwith trial nets"
                                 ), 
                      values = c("black", "darkgreen", "skyblue"),
                      name = "",
                      breaks = c("None", "IG2", "IG2_tu")) +
  scale_fill_manual(labels = c("No trial nets added",
                                 "Trial nets replaced\nwith pyrethroid-only nets",
                                 "Trial nets replaced\nwith trial nets"
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
              alpha = 0.4) +
  
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
                                 "Trial nets replaced\nwith trial nets"
  ), 
  values = c("black", "darkgreen", "skyblue"),
  name = "",
  breaks = c("None", "IG2", "IG2_tu")) +
  scale_fill_manual(labels = c("No trial nets added",
                               "Trial nets replaced\nwith pyrethroid-only nets",
                               "Trial nets replaced\nwith trial nets"
  ), 
  values = c("black", "darkgreen", "skyblue"),
  name = "",
  breaks = c("None", "IG2", "IG2_tu"))

##### actual vs fitted plot #####
af_arm_sims_T_tu <- data.frame("t" = rep(pred_prev_T[[1]]$timestep, 3),
                            "date" = rep(pred_prev_T[[1]]$date, 3),
                            "net" = c(rep("IG2", sim_length),
                                      rep("PBO", sim_length),
                                      rep("Pyrethroid_only", sim_length)),
                            "prev" = c(extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                                                   Country_in = "Tanzania", bioassay_uncertainty = "middle",
                                                                   vals_pred = T_vals_pred_tu, pred_prev = pred_prev_T_tu),
                                                      
                                                      extract_prev(Net_in = "PBO", int_Net_in = "PBO",
                                                                   Country_in = "Tanzania", bioassay_uncertainty = "middle",
                                                                   vals_pred = T_vals_pred_tu, pred_prev = pred_prev_T_tu),
                                                      
                                                      extract_prev(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only",
                                                                   Country_in = "Tanzania", bioassay_uncertainty = "middle",
                                                                   vals_pred = T_vals_pred_tu, pred_prev = pred_prev_T_tu)
                                        ),
                                        
                                        "lower" = c(extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                                                 Country_in = "Tanzania", bioassay_uncertainty = "lower",
                                                                 vals_pred = T_vals_pred_tu, pred_prev = pred_prev_T_tu),
                                                    
                                                    extract_prev(Net_in = "PBO", int_Net_in = "PBO",
                                                                 Country_in = "Tanzania", bioassay_uncertainty = "lower",
                                                                 vals_pred = T_vals_pred_tu, pred_prev = pred_prev_T_tu),
                                                    
                                                    extract_prev(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only",
                                                                 Country_in = "Tanzania", bioassay_uncertainty = "lower",
                                                                 vals_pred = T_vals_pred_tu, pred_prev = pred_prev_T_tu)
                                        ),
                                        
                                        "upper" = c(extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                                                 Country_in = "Tanzania", bioassay_uncertainty = "upper",
                                                                 vals_pred = T_vals_pred_tu, pred_prev = pred_prev_T_tu),
                                                    
                                                    extract_prev(Net_in = "PBO", int_Net_in = "PBO",
                                                                 Country_in = "Tanzania", bioassay_uncertainty = "upper",
                                                                 vals_pred = T_vals_pred_tu, pred_prev = pred_prev_T_tu),
                                                    
                                                    extract_prev(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only",
                                                                 Country_in = "Tanzania", bioassay_uncertainty = "upper",
                                                                 vals_pred = T_vals_pred_tu, pred_prev = pred_prev_T_tu)
                                        )
)
                                        
af_arm_sims_B_tu <- data.frame("t" = rep(pred_prev_B[[1]]$timestep, 2),
                            "date" = rep(pred_prev_B[[1]]$date, 2),
                            "net" = c(rep("IG2", sim_length),
                                      rep("Pyrethroid_only", sim_length)),
                            "prev" = c(extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                                    Country_in = "Benin", bioassay_uncertainty = "middle",
                                                    vals_pred = B_vals_pred_tu, pred_prev = pred_prev_B_tu),
                                       
                                       extract_prev(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only",
                                                    Country_in = "Benin", bioassay_uncertainty = "middle",
                                                    vals_pred = B_vals_pred_tu, pred_prev = pred_prev_B_tu)
                            ),
                            "lower" = c(extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                                     Country_in = "Benin", bioassay_uncertainty = "lower",
                                                     vals_pred = B_vals_pred_tu, pred_prev = pred_prev_B_tu),
                                        
                                        extract_prev(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only",
                                                     Country_in = "Benin", bioassay_uncertainty = "lower",
                                                     vals_pred = B_vals_pred_tu, pred_prev = pred_prev_B_tu)
                            ),
                            "upper" = c(extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                                     Country_in = "Benin", bioassay_uncertainty = "upper",
                                                     vals_pred = B_vals_pred_tu, pred_prev = pred_prev_B_tu),
                                        
                                        extract_prev(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only",
                                                     Country_in = "Benin", bioassay_uncertainty = "upper",
                                                     vals_pred = B_vals_pred_tu, pred_prev = pred_prev_B_tu)
                            )
)

inds_T <- match(interaction(df_T$date, df_T$Net), interaction(arm_sims_T$date, arm_sims_T$net))

inds_T_tu <- match(interaction(df_T$date, df_T$Net), interaction(af_arm_sims_T_tu$date, af_arm_sims_T_tu$net))

df_T <- as.data.frame(df_T) %>% mutate(pred_prev = arm_sims_T[inds_T, "prev"],
                                       pred_prev_l = arm_sims_T[inds_T, "lower"],
                                       pred_prev_u = arm_sims_T[inds_T, "upper"],
                                       pred_prev_cf = af_arm_sims_T_tu[inds_T_tu, "prev"],
                                       pred_prev_l_cf = af_arm_sims_T_tu[inds_T_tu, "lower"],
                                       pred_prev_u_cf = af_arm_sims_T_tu[inds_T_tu, "upper"])

inds_B <- match(interaction(df_B$date, df_B$Net), interaction(arm_sims_B$date, arm_sims_B$net))

inds_B_tu <- match(interaction(df_B$date, df_B$Net), interaction(af_arm_sims_B_tu$date, af_arm_sims_B_tu$net))

df_B <- as.data.frame(df_B) %>% mutate(pred_prev = arm_sims_B[inds_B, "prev"],
                                       pred_prev_l = arm_sims_B[inds_B, "lower"],
                                       pred_prev_u = arm_sims_B[inds_B, "upper"],
                                       pred_prev_cf = af_arm_sims_B_tu[inds_B_tu, "prev"],
                                       pred_prev_l_cf = af_arm_sims_B_tu[inds_B_tu, "lower"],
                                       pred_prev_u_cf = af_arm_sims_B_tu[inds_B_tu, "upper"])

df_af <- rbind(df_T[, c("Net", "Country", "date", "Malaria_prevalence", "Malaria_prevalence_u", 
                        "Malaria_prevalence_l", "pred_prev", "pred_prev_l", "pred_prev_u")], 
               df_B[, c("Net", "Country", "date", "Malaria_prevalence", "Malaria_prevalence_u",
                        "Malaria_prevalence_l", "pred_prev", "pred_prev_l", "pred_prev_u")])

df_af <- df_af[is.na(df_af$Malaria_prevalence)!=1, ]

af_plot <- ggplot(data = subset(df_af, !(date %in% c(baseline_start_date_B,
                                                     baseline_start_date_T)) & Net != "RG")) +
  geom_abline(slope = 1, linetype = 2, linewidth = 1.1) +
  # geom_smooth(formula = y ~ x-1, se = TRUE, aes(x = Malaria_prevalence, y = pred_prev),
  #              method = "lm", col = "grey50", fullrange = TRUE, alpha = 0.1) +
  geom_pointrange(size = 1, 
                  aes(x = Malaria_prevalence, y = pred_prev, ymin = pred_prev_u, ymax = pred_prev_l, col = Net, shape = Country)) +
  geom_errorbarh(aes(y = pred_prev, xmin = Malaria_prevalence_l, xmax = Malaria_prevalence_u, col = Net, group = Country)) +
  theme_classic() + 
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 14),
        legend.title = element_text(size=14),
        legend.background = element_blank(),
        legend.position = c(0.85, 0.2),
        legend.box.background = element_rect(colour = "black")) +
  scale_colour_manual(values = c("blue", "aquamarine", "darkgreen"),
                      labels = c("pyrethroid-only", "pyrethroid-PBO", "pyrethroid-pyrrole"),
                      breaks =c("Pyrethroid_only", "PBO", "IG2")) +
  scale_shape_manual(values = c(15, 16)) +
  ylab("Predicted malaria prevalence") + xlab("Observed malaria prevalence") +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent, breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent, breaks = seq(0, 1, 0.1)) +
  coord_cartesian(xlim = c(0, 0.7), ylim = c(0, 0.7))




##### extracting the results #####

calc_incidence <- function(output, s_date, e_date, min_age, max_age){
  s_time <- which(output$date == s_date)
  e_time <- which(output$date == e_date)
  return(output[s_time:e_time,paste0("n_inc_clinical_",min_age,"_",max_age)]/
           output[s_time:e_time, paste0("n_", min_age,"_", max_age)])
}

extract_inc <- function(Net_in, 
                        int_Net_in,
                        bioassay_uncertainty,
                        vals_pred, 
                        pred_prev,
                        top_up_net = NA,
                        n_rep = 3,
                        s_date = intervention_date_T, 
                        e_date = intervention_date_T + 365,
                        min_age = 0.5*365,
                        max_age = 10*365){
  
  #p <- if(bioassay_uncertainty == "middle"){0.5} else if(bioassay_uncertainty == "lower"){0.025} else{0.975}
  
  inds <- if(is.na(top_up_net) == 1){
    which(vals_pred$Net_in == Net_in & 
            vals_pred$int_Net_in == int_Net_in &
            vals_pred$bioassay_uncertainty == bioassay_uncertainty)
  } else{
    which(vals_pred$Net_in == Net_in & 
            vals_pred$int_Net_in == int_Net_in &
            vals_pred$bioassay_uncertainty == bioassay_uncertainty &
            vals_pred$top_up_net == top_up_net)
  }
  
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
    
    return(mean(mat)) # incidence is per person
}

# calculate efficacy - put on a different plot
# put grey points from previous paper - R2 with and without the new ones
# R2 from the x=y line
rel_reduction <- rbind(data.frame(s_times = rep(c(intervention_date_T, 
                                                  intervention_date_T + 365,
                                                  intervention_date_T,
                                                  intervention_date_T),
                                                3),
                            e_times = rep(c(intervention_date_T + 365,
                                            intervention_date_T + 365 * 2,
                                            intervention_date_T + 365 * 2,
                                            intervention_date_T + 365 * 3),
                                          3),
                            year = rep(c("Year 1", "Year 2", "Overall", "0_to_3"), 3),
                            net = sort(rep(c("Pyrethroid_only", "PBO", "IG2"), 4))
) %>% rowwise() %>% 
  mutate(p_only_m = extract_inc(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only",
                                bioassay_uncertainty = "middle",
                                vals_pred = T_vals_pred,
                                pred_prev = pred_prev_T,
                                s_date = s_times,
                                e_date = e_times),
         p_only_l = extract_inc(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only",
                                bioassay_uncertainty = "lower",
                                vals_pred = T_vals_pred,
                                pred_prev = pred_prev_T,
                                s_date = s_times,
                                e_date = e_times),
         p_only_u = extract_inc(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only",
                                bioassay_uncertainty = "upper",
                                vals_pred = T_vals_pred,
                                pred_prev = pred_prev_T,
                                s_date = s_times,
                                e_date = e_times),
         t_net_m = extract_inc(Net_in = net, int_Net_in = net,
                      bioassay_uncertainty = "middle",
                      vals_pred = T_vals_pred,
                      pred_prev = pred_prev_T,
                      s_date = s_times,
                      e_date = e_times),
         t_net_l = extract_inc(Net_in = net, int_Net_in = net,
                               bioassay_uncertainty = "lower",
                               vals_pred = T_vals_pred,
                               pred_prev = pred_prev_T,
                               s_date = s_times,
                               e_date = e_times),
         t_net_u = extract_inc(Net_in = net, int_Net_in = net,
                               bioassay_uncertainty = "upper",
                               vals_pred = T_vals_pred,
                               pred_prev = pred_prev_T,
                               s_date = s_times,
                               e_date = e_times),
         
         tu_net_m = extract_inc(Net_in = net, int_Net_in = net,
                                bioassay_uncertainty = "middle",
                                vals_pred = T_vals_pred_tu,
                                pred_prev = pred_prev_T_tu,
                                s_date = s_times,
                                e_date = e_times),
         
         tu_net_l = extract_inc(Net_in = net, int_Net_in = net,
                                bioassay_uncertainty = "lower",
                                vals_pred = T_vals_pred_tu,
                                pred_prev = pred_prev_T_tu,
                                s_date = s_times,
                                e_date = e_times),
         
         tu_net_u = extract_inc(Net_in = net, int_Net_in = net,
                                bioassay_uncertainty = "upper",
                                vals_pred = T_vals_pred_tu,
                                pred_prev = pred_prev_T_tu,
                                s_date = s_times,
                                e_date = e_times),
         
         cf_net_m = extract_inc(Net_in = net, int_Net_in = net,
                               bioassay_uncertainty = "middle",
                               vals_pred = T_vals_pred_tu_none,
                               pred_prev = pred_prev_T_tu_none,
                               s_date = s_times,
                               e_date = e_times),
         
         cf_net_l = extract_inc(Net_in = net, int_Net_in = net,
                               bioassay_uncertainty = "lower",
                               vals_pred = T_vals_pred_tu_none,
                               pred_prev = pred_prev_T_tu_none,
                               s_date = s_times,
                               e_date = e_times),
         
         cf_net_u = extract_inc(Net_in = net, int_Net_in = net,
                               bioassay_uncertainty = "upper",
                               vals_pred = T_vals_pred_tu_none,
                               pred_prev = pred_prev_T_tu_none,
                               s_date = s_times,
                               e_date = e_times),
         
    Country = "Tanzania"
  ),
data.frame(s_times = rep(c(intervention_date_B, 
                           intervention_date_B + 365,
                           intervention_date_B,
                           intervention_date_B), 2),
           e_times = rep(c(intervention_date_B + 365,
                           intervention_date_B + 365 * 2,
                           intervention_date_B + 365 * 2,
                           intervention_date_B + 365 * 3), 2),
           year = rep(c("Year 1", "Year 2", "Overall", "0_to_3"), 2),
           net = sort(rep(c("Pyrethroid_only", "IG2"), 4))) %>% 
  rowwise() %>% 
  mutate(p_only_m = extract_inc(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only",
                                bioassay_uncertainty = "middle",
                                vals_pred = B_vals_pred,
                                pred_prev = pred_prev_B,
                                s_date = s_times,
                                e_date = e_times),
         p_only_l = extract_inc(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only",
                                bioassay_uncertainty = "lower",
                                vals_pred = B_vals_pred,
                                pred_prev = pred_prev_B,
                                s_date = s_times,
                                e_date = e_times),
         p_only_u = extract_inc(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only",
                                bioassay_uncertainty = "upper",
                                vals_pred = B_vals_pred,
                                pred_prev = pred_prev_B,
                                s_date = s_times,
                                e_date = e_times),
         t_net_m = extract_inc(Net_in = net, int_Net_in = net,
                               bioassay_uncertainty = "middle",
                               vals_pred = B_vals_pred,
                               pred_prev = pred_prev_B,
                               s_date = s_times,
                               e_date = e_times),
         t_net_l = extract_inc(Net_in = net, int_Net_in = net,
                               bioassay_uncertainty = "lower",
                               vals_pred = B_vals_pred,
                               pred_prev = pred_prev_B,
                               s_date = s_times,
                               e_date = e_times),
         t_net_u = extract_inc(Net_in = net, int_Net_in = net,
                               bioassay_uncertainty = "upper",
                               vals_pred = B_vals_pred,
                               pred_prev = pred_prev_B,
                               s_date = s_times,
                               e_date = e_times),
         
         tu_net_m = extract_inc(Net_in = net, int_Net_in = net,
                     bioassay_uncertainty = "middle",
                     vals_pred = B_vals_pred_tu,
                     pred_prev = pred_prev_B_tu,
                     s_date = s_times,
                     e_date = e_times),
         
         tu_net_l = extract_inc(Net_in = net, int_Net_in = net,
                                bioassay_uncertainty = "lower",
                                vals_pred = B_vals_pred_tu,
                                pred_prev = pred_prev_B_tu,
                                s_date = s_times,
                                e_date = e_times),
         
         tu_net_u = extract_inc(Net_in = net, int_Net_in = net,
                                bioassay_uncertainty = "upper",
                                vals_pred = B_vals_pred_tu,
                                pred_prev = pred_prev_B_tu,
                                s_date = s_times,
                                e_date = e_times),
         
         cf_net_m = extract_inc(Net_in = net, int_Net_in = net,
                                bioassay_uncertainty = "middle",
                                vals_pred = B_vals_pred_tu_none,
                                pred_prev = pred_prev_B_tu_none,
                                s_date = s_times,
                                e_date = e_times),
         
         cf_net_l = extract_inc(Net_in = net, int_Net_in = net,
                                bioassay_uncertainty = "lower",
                                vals_pred = B_vals_pred_tu_none,
                                pred_prev = pred_prev_B_tu_none,
                                s_date = s_times,
                                e_date = e_times),
         
         cf_net_u = extract_inc(Net_in = net, int_Net_in = net,
                                bioassay_uncertainty = "upper",
                                vals_pred = B_vals_pred_tu_none,
                                pred_prev = pred_prev_B_tu_none,
                                s_date = s_times,
                                e_date = e_times),
         
         Country = "Benin")
) %>% rowwise() %>% mutate(rr_m = (p_only_m - t_net_m)/p_only_m,
                           rr_l = (p_only_l - t_net_l)/p_only_l,
                           rr_u = (p_only_u - t_net_u)/p_only_u,
                           
                           overall_efficacy_m = (cf_net_m - t_net_m) / cf_net_m,
                           overall_efficacy_l = (cf_net_l - t_net_l) / cf_net_l,
                           overall_efficacy_u = (cf_net_u - t_net_u) / cf_net_u,
                           
                           rr_tu_m = (p_only_m - tu_net_m)/p_only_m,
                           rr_tu_l = (p_only_l - tu_net_l)/p_only_l,
                           rr_tu_u = (p_only_u - tu_net_u)/p_only_u,
                           
                           overall_efficacy_tu_m = (cf_net_m - tu_net_m) / cf_net_m,
                           overall_efficacy_tu_l = (cf_net_l - tu_net_l) / cf_net_l,
                           overall_efficacy_tu_u = (cf_net_u - tu_net_u) / cf_net_u
                           
                           )

actual_inc <- read_xlsx("data/actual_incidence_estimates.xlsx") %>% mutate(rr_m = (Pyrethroid_only_incidence - Trial_incidence)/Pyrethroid_only_incidence,
                                                                           sample = "observed")

rel_reduction$Country <- factor(rel_reduction$Country, levels=c("Tanzania", "Benin"))
actual_inc$Country <- factor(actual_inc$Country, levels=c("Tanzania", "Benin"))

rr_plot <- ggplot(data = rbind(subset(rel_reduction, net != "Pyrethroid_only" & year != "0_to_3")[,c("Country", "year", "rr_m", "net")] %>% 
                                 mutate(sample = "predicted"),
                    actual_inc[,c("Country", "year", "rr_m", "net", "sample")]) %>% rowwise() %>%  
                    mutate(net = if(net == "IG2"){"pyrethroid-pyrrole"} else{"pyrethroid-PBO"}), 
       aes(x = year, y = rr_m, fill = net, group = interaction(net, sample))) + 
  ggpattern::geom_bar_pattern(stat = "identity", 
                   alpha = 0.4, 
                   position = position_dodge(preserve = "single"),
                   pattern_colour = "black",
                   col = "black",
                   aes(pattern = sample)) +
  scale_pattern_manual(values = c(predicted = "circle", observed = "none"), name = "Sample") +
  facet_wrap(~Country + net) +
  xlab("Year post trial net distribution") + ylab("Efficacy (clinical incidence averted relative\nto the value in the pyrethroid-only arm)") +
  theme_classic() + theme(text = element_text(size = 20),
                          axis.text.y = element_text(size = 17),
                          axis.text.x = element_text(size = 12.75),
                          legend.text = element_text(size = 14),
                          legend.title = element_text(size=14),
                          legend.background = element_blank()) +
  scale_y_continuous(labels = scales::percent) +
  scale_colour_manual(values = c("aquamarine", "darkgreen"), name = "Net") +
  scale_fill_manual(values = c("aquamarine", "darkgreen"), name = "Net")

png(file = "figures/model_simulations.png", height = 1100, width = 1800)

(T_arm_plot | B_arm_plot | af_plot) / (rr_plot | cf_plot_T | cf_plot_B) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold'))
    
dev.off()

############################################
##### extracting the prevalence values #####
############################################

# extracting the values
inc_out <- rel_reduction %>% rowwise() %>%
  mutate(observed_arm = paste0(round(t_net_m, digits = 2), " (", round(rr_m, digits = 2)*100,"%)"),
         overall_arm = paste0(round(overall_efficacy_m, digits = 2)*100, "%"),
         observed_tu = paste0(round(tu_net_m, digits = 2), " (", round(rr_tu_m, digits = 2)*100,"%)"),
         overall_tu = paste0(round(overall_efficacy_tu_m, digits = 2)*100, "%")) %>% 
  select(c(Country, year, net, observed_arm, overall_arm, observed_tu, overall_tu))

write.csv(inc_out, file = "data/inc_estimates.csv")

write.csv(actual_inc %>% mutate(inc_out = paste0(Trial_incidence, " (", round(rr_m, digits = 2)*100, "%)")), file = "data/actual_inc_estimates.csv")

# point prevalence values at the given times

df_T <- left_join(df_T, data.frame(date = rep(pred_prev_T_tu_none[[1]]$date, 3),
                           Net = sort(rep(c("IG2", "PBO", "Pyrethroid_only"), length(pred_prev_T_tu_none[[1]]$date))),
                           prev_none = c(
                             extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                          Country_in = "Tanzania", bioassay_uncertainty = "middle",
                                          vals_pred = T_vals_pred_tu_none, pred_prev = pred_prev_T_tu_none),
                             
                             extract_prev(Net_in = "PBO", int_Net_in = "PBO",
                                          Country_in = "Tanzania", bioassay_uncertainty = "middle",
                                          vals_pred = T_vals_pred_tu_none, pred_prev = pred_prev_T_tu_none),
                             
                             extract_prev(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only",
                                          Country_in = "Tanzania", bioassay_uncertainty = "middle",
                                          vals_pred = T_vals_pred_tu_none, pred_prev = pred_prev_T_tu_none))), 
          by = c("date", "Net"))

df_B <- left_join(df_B, data.frame(date = rep(pred_prev_B_tu_none[[1]]$date, 2),
                                   Net = sort(rep(c("IG2", "Pyrethroid_only"), length(pred_prev_B_tu_none[[1]]$date))),
                                   prev_none = c(
                                     extract_prev(Net_in = "IG2", int_Net_in = "IG2",
                                                  Country_in = "Benin", bioassay_uncertainty = "middle",
                                                  vals_pred = B_vals_pred_tu_none, pred_prev = pred_prev_B_tu_none),
                                     
                                     extract_prev(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only",
                                                  Country_in = "Benin", bioassay_uncertainty = "middle",
                                                  vals_pred = B_vals_pred_tu_none, pred_prev = pred_prev_B_tu_none))), 
                  by = c("date", "Net"))


df_out <- rbind(left_join(subset(df_T, Time_months %in% c(6, 12, 18, 24) & Net!= "RG")[, c("Net", "Time_months", "Malaria_prevalence", "pred_prev",
                                                                    "pred_prev_cf", "prev_none")],
                          
                          subset(df_T, Time_months %in% c(6, 12, 18, 24) & Net == "Pyrethroid_only")[, c("Time_months", "Malaria_prevalence",  "pred_prev",
                                                                           "pred_prev_cf", "prev_none")] %>% 
                            setNames(c("Time_months", "Malaria_prevalence_p_only", "pred_prev_p_only", "pred_prev_cf_p_only", "prev_none_p_only")),
          by = c("Time_months")) %>% 
            arrange(Net, Time_months) %>% mutate(Country = "Tanzania"),
  
  left_join(subset(df_B, Time_months %in% c(6, 12, 18, 24) & Net!= "RG")[, c("Net", "Time_months", "Malaria_prevalence", "pred_prev",
                                                                            "pred_prev_cf", "prev_none")],
            
           subset(df_B, Time_months %in% c(6, 12, 18, 24) & Net == "Pyrethroid_only")[, c("Time_months", "Malaria_prevalence", "pred_prev", "pred_prev_cf", "prev_none")] %>% 
             setNames(c("Time_months", "Malaria_prevalence_p_only", "pred_prev_p_only", "pred_prev_cf_p_only", "prev_none_p_only")),
           by = c("Time_months")) %>% 
                  arrange(Net, Time_months) %>% mutate(Country = "Benin")
  ) %>% 
  mutate(observed = paste0(round(Malaria_prevalence, digits = 2), 
                                    " (", round((Malaria_prevalence_p_only - Malaria_prevalence) / Malaria_prevalence_p_only, digits = 2)*100,"%)"),
         
         model_observed = paste0(round(pred_prev, digits = 2), 
                                       " (", round((pred_prev_p_only - pred_prev)/pred_prev_p_only, digits = 2)*100,"%)"),
         
         overall = paste0(round((prev_none - pred_prev)/prev_none, digits = 2)*100, "%"),
         
         model_observed_tu = paste0(round(pred_prev_cf, digits = 2), 
                                 " (", round((pred_prev_p_only - pred_prev_cf)/pred_prev_p_only, digits = 2)*100,"%)"),
         
         overall_tu = paste0(round((prev_none - pred_prev_cf)/prev_none, digits = 2)*100, "%")
         )

write.csv(df_out, file = "data/prevalence_estimates.csv")

# annual mean prevalence values

extract_prev_mean_annual <- function(Net_in, 
                         int_Net_in,
                         Country_in,
                         bioassay_uncertainty = "middle",
                         vals_pred, 
                         pred_prev,
                         top_up_net = NA,
                         n_rep = 3,
                         trial_net_cov_in = NULL,
                         s_date,
                         e_date){
  
  #p <- if(bioassay_uncertainty == "middle"){0.5} else if(bioassay_uncertainty == "lower"){0.025} else{0.975}
  
  if(is.null(trial_net_cov_in) == 1){
    inds <- if(is.na(top_up_net) == 1){
      which(vals_pred$Net_in == Net_in & 
              vals_pred$int_Net_in == int_Net_in &
              vals_pred$bioassay_uncertainty == bioassay_uncertainty)
    } else{
      which(vals_pred$Net_in == Net_in & 
              vals_pred$int_Net_in == int_Net_in &
              vals_pred$bioassay_uncertainty == bioassay_uncertainty &
              vals_pred$top_up_net == top_up_net)
    }} else{
      
      inds <- if(is.na(top_up_net) == 1){
        which(vals_pred$Net_in == Net_in & 
                vals_pred$int_Net_in == int_Net_in &
                vals_pred$bioassay_uncertainty == bioassay_uncertainty &
                is.na(vals_pred$trial_net_cov_in) == is.na(trial_net_cov_in))
      } else{
        which(vals_pred$Net_in == Net_in & 
                vals_pred$int_Net_in == int_Net_in &
                vals_pred$bioassay_uncertainty == bioassay_uncertainty &
                vals_pred$top_up_net == top_up_net &
                is.na(vals_pred$trial_net_cov_in) == is.na(trial_net_cov_in))
      }
      
      
    }
  
  if(length(inds) != n_rep){return(NA)}else{
    
    prev_mat <- if(Country_in == "Tanzania"){sapply(inds, function(i, p){summary_pfpr_0.5_5(p[[i]])}, p = pred_prev)
    } else if(Country_in == "Benin"){sapply(inds, function(i, p){summary_pfpr_0_100_all(p[[i]])}, p = pred_prev)
    }
    
    s_time <- which(pred_prev[[1]]$date == s_date)
    e_time <- which(pred_prev[[1]]$date == e_date)
    
    
    prev_mat <- apply(prev_mat[s_time:e_time,], 2, mean)
    
    return(mean(prev_mat))
  }
}

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

