rm(list = ls())

# find proportion of people with bednets

# turning on the development mode for the model with LSM

suppressPackageStartupMessages(library(ggplot2)); library(malariasimulation); library(malariaEquilibrium)
library(reshape2); library(tidyverse); library(readxl); library(rstan); library(pracma)
library(foresite); library(doParallel); library(foreach);
library(cowplot); library(pammtools);

source(file = "functions.R"); source(file = "retention_fit_top_up.R")

### functions
# can drop the coverage by each time-frame
x <- decline_ds(seq(0, 100), 0.5, 0.01)

for(i in 2:length(x)){
  print(decline_ds(1, dn0 = x[i-1], gamma_n = 0.01))
}

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
# include a warm up period of ~40 years
sim_length <- 50 * 365

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

df_B <- rbind(subset(df, Country == "Benin" & Time_months == 0) %>% mutate(date = baseline_start_date_B),
              subset(df, Country == "Benin" & Time_months != 0) %>% mutate(date = as.Date(Time_months*30 + intervention_date_B)))

df_T <- rbind(subset(df, Country == "Tanzania" & Time_months == 0) %>% mutate(date = baseline_start_date_T),
              subset(df, Country == "Tanzania" & Time_months != 0) %>% mutate(date = as.Date(Time_months*30 + intervention_date_T)))

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
                          tu_diff_time = 1,
                          sim_length = sim_length,
                          int_bed_net_time = int_time_B)

top_up_p_only_B <- get_top_up(net = "Pyrethroid_only",
                            df = df,
                            country = "Benin",
                            tu_diff_time = 1,
                            sim_length = sim_length,
                            int_bed_net_time = int_time_B)

top_up_IG2_T_0.25y <- get_top_up(net = "IG2",
                           df = df,
                           country = "Tanzania",
                           tu_diff_time = 364/4,
                           sim_length = sim_length,
                           int_bed_net_time = int_time_T)

top_up_IG2_T_1y <- get_top_up(net = "IG2",
                           df = df,
                           country = "Tanzania",
                           tu_diff_time = 364,
                           sim_length = sim_length,
                           int_bed_net_time = int_time_T)

top_up_IG2_T_3y <- get_top_up(net = "IG2",
                              df = df,
                              country = "Tanzania",
                              tu_diff_time = 365*3,
                              sim_length = sim_length,
                              int_bed_net_time = int_time_T)


top_up_RG_T <- get_top_up(net = "RG",
                           df = df,
                           country = "Tanzania",
                           tu_diff_time = 1,
                           sim_length = sim_length,
                           int_bed_net_time = int_time_T)

top_up_PBO_T <- get_top_up(net = "PBO",
                           df = df,
                           country = "Tanzania",
                           tu_diff_time = 1,
                           sim_length = sim_length,
                           int_bed_net_time = int_time_T)

top_up_p_only_T <- get_top_up(net = "Pyrethroid_only",
                           df = df,
                           country = "Tanzania",
                           tu_diff_time = 1,
                           sim_length = sim_length,
                           int_bed_net_time = int_time_T)

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
png(file = "figures/net_cover.png", height = 750, width = 1100)
plot_grid(
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
      
      #model_seasonality = TRUE, ## Seasonality to match study site inputs [sites_13]
      model_seasonality = FALSE,
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
  
  init_cov <- min(df_in[baseline_index, "Bed_net_use_both"], top_up$covs_df[1,"mean_cov_l"]) #top_up$mean_bed_net_use_both # changed to baseline coverage
  
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
clusterExport(cl = cl, varlist = c("T_site", "int_time_T", "baseline_time_T", "df", "df_species_in", 
                                   "dat_res_pyr", "sim_length", 
                                   "top_up_IG2_T", "top_up_PBO_T", 
                                   "top_up_RG_T", "top_up_p_only_T"))

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
    sim_years = sim_length/365)},
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
clusterExport(cl = cl, varlist = c("B_site", "int_time_B", "baseline_time_B", "df", "df_species_in", 
                                   "dat_res_pyr", "sim_length", 
                                   "top_up_IG2_B",
                                   "top_up_RG_B", "top_up_p_only_B"))
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
    sim_years = sim_length/365)},
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
                        bioassay_uncertainty = "middle"){
  
  top_up <- get(top_up_name)
  
  # specifying the retention parameter
  # if no top up then let the nets decline
  retention <- if(top_up_net == "none"){mean(rstan::extract(top_up$fit_base, "retention")$retention) * 365}else{10^10/365}
  
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
  gamman <- dat_res_pyr_in[bioassay_index, gamman_name]
  
  dat_int <- if(int_Net_in == "IG2"){dat_res_pp_in} else if(int_Net_in == "PBO"){dat_res_pbo_in} else if(int_Net_in == "PPF"){dat_res_ppf_in} else if(int_Net_in == "Pyrethroid_only"){dat_res_pyr_in}
    
  int_bioassay_index <- which(round(dat_int$bioassay_mortality, digits = 2) == bioassay_mortality)
  dn0_int <- dat_int[int_bioassay_index, dn0_name]
  rn0_int <- dat_int[int_bioassay_index, rn0_name]
  gamman_int <- dat_int[int_bioassay_index, gamman_name]
  
  dat_top_up <- if(top_up_net == "IG2"){dat_res_pp_in}else if(top_up_net == "PBO"){dat_res_pbo_in} else if(top_up_net == "PPF"){dat_res_ppf_in} else if(top_up_net == "Pyrethroid_only"){dat_res_pyr_in}
  
  if(top_up_net != "none"){
    tu_bioassay_index <- which(round(dat_top_up$bioassay_mortality, digits = 2) == bioassay_mortality)
    dn0_tu <- dat_top_up[tu_bioassay_index, dn0_name] #if(top_up_net == "IG2"){0.85}else if(top_up_net == "Pyrethroid_only"){0.05}
    rn0_tu <- dat_top_up[tu_bioassay_index, rn0_name] #if(top_up_net == "IG2"){0.13}else if(top_up_net == "Pyrethroid_only"){0.05}
    gamman_tu <- dat_top_up[tu_bioassay_index, gamman_name]
  }
  
  
  bednetparams <- simparams
  
  init_cov <- min(df_in[baseline_index, "Bed_net_use_both"], top_up$covs_df[1,"mean_cov_l"]) # top_up$mean_bed_net_use_both
  covs_df <- top_up$covs_df
  
  if(top_up_net == "none"){
    init_times <- round(seq(1, intervention_time, 365*5), digits = 0)
    n_init_times <- length(init_times)
    
    bednet_events = data.frame(timestep = c(init_times, intervention_time),
                               name=c(paste0(rep("baseline_nets", n_init_times),
                                             ": ", init_times), 
                                      "intervention_nets"))
    
    bednetparams_1 <- set_bednets(bednetparams,
                                  
                                  timesteps = bednet_events$timestep,
                                  
                                  coverages = c(rep(top_up$mean_bed_net_use_both, n_init_times), top_up$mean_bed_net_use_both),
                                  
                                  retention = retention, # assumed
                                  
                                  # each row needs to show the efficacy parameter across years (and cols are diff mosquito)
                                  # gambiae, coluzzi, arabiensis, funestus, coustani
                                  # no resistance assumed for coustani
                                  
                                  dn0 = matrix(c(rep(dn0, n_init_times), dn0_int, 
                                                 rep(dn0, n_init_times), dn0_int, 
                                                 rep(dn0, n_init_times), dn0_int), nrow = n_init_times + 1, ncol = 3),
                                  rn = matrix(c(rep(rn0, n_init_times), rn0_int, 
                                                rep(rn0, n_init_times), rn0_int, 
                                                rep(rn0, n_init_times), rn0_int), nrow = n_init_times + 1, ncol = 3),
                                  rnm = matrix(c(rep(rnm, n_init_times), rnm,
                                                 rep(rnm, n_init_times), rnm, 
                                                 rep(rnm, n_init_times), rnm), nrow = n_init_times + 1, ncol = 3),
                                  gamman = c(rep(gamman, n_init_times), gamman_int) * 365
    )
    
  } else{
    
  n_top_up <- nrow(covs_df) - 1
  
  # less top up if bioassay mortality is higher
  top_up_name <- case_when(bioassay_uncertainty == "middle" ~ "top_up",
                           bioassay_uncertainty == "lower" ~ "top_up_l",
                           bioassay_uncertainty == "upper" ~ "top_up_u"
                           )
  
  init_cov_name <- case_when(bioassay_uncertainty == "middle" ~ "mean_cov_m",
                             bioassay_uncertainty == "lower" ~ "mean_cov_l",
                             bioassay_uncertainty == "upper" ~ "mean_cov_u")
  
  # initial bednet event is at time 1
  bednet_events = data.frame(timestep = c(1, covs_df$t + intervention_time),
                             name=c("baseline_nets", "intervention_nets", paste0("top_up_nets: ", seq(1, nrow(covs_df)-1)))
  )
  
  bednetparams_1 <- set_bednets(bednetparams,

                                timesteps = bednet_events$timestep,

                                coverages = c(init_cov, covs_df[1, init_cov_name], covs_df[-1, top_up_name]), # initial coverage is the mean coverage between first two

                                retention = retention, # assumed

                                # each row needs to show the efficacy parameter across years (and cols are diff mosquito)
                                # gambiae, coluzzi, arabiensis, funestus, coustani
                                # no resistance assumed for coustani

                                dn0 = matrix(c(dn0, dn0_int, rep(dn0_tu, n_top_up), 
                                               dn0, dn0_int, rep(dn0_tu, n_top_up),
                                               dn0, dn0_int, rep(dn0_tu, n_top_up)), nrow = n_top_up + 2, ncol = 3),
                                rn = matrix(c(rn0, rn0_int, rep(rn0_tu, n_top_up), 
                                              rn0, rn0_int, rep(rn0_tu, n_top_up),
                                              rn0, rn0_int, rep(rn0_tu, n_top_up)), nrow = n_top_up + 2, ncol = 3),
                                rnm = matrix(c(rnm, rnm, rep(rnm, n_top_up), 
                                               rnm, rnm, rep(rnm, n_top_up),
                                               rnm, rnm, rep(rnm, n_top_up)), nrow = n_top_up + 2, ncol = 3),
                                gamman = c(gamman * 365, 1000000, rep(1000000, n_top_up))#c(gamman, gamman_int, rep(gamman_tu, n_top_up)) * 365
  )
    
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

out_none <- sim_forward(start_EIR = T_vals[1, "start_EIR"],
            
            Location_in = T_vals[1, "Location_in"],
            Net_in = T_vals[1, "Net_in"],
            Country_in = T_vals[1, "Country_in"],
            
            int_Net_in = T_vals[1, "Net_in"],
            
            top_up_name = T_vals[1, "top_up_name"],
            
            top_up_net = "none",
            
            sites_in = T_site,
            
            intervention_time = int_time_T,
            baseline_pred_days = baseline_time_T)

out_p_0.25 <- sim_forward(start_EIR = T_vals[1, "start_EIR"],
                        
                        Location_in = T_vals[1, "Location_in"],
                        Net_in = T_vals[1, "Net_in"],
                        Country_in = T_vals[1, "Country_in"],
                        
                        int_Net_in = T_vals[1, "Net_in"],
                        
                        top_up_name = "top_up_IG2_T_0.25y",
                        
                        top_up_net = "Pyrethroid_only",
                        
                        sites_in = T_site,
                        
                        intervention_time = int_time_T,
                        baseline_pred_days = baseline_time_T)

out_IG2_0.25 <- sim_forward(start_EIR = T_vals[1, "start_EIR"],
                       
                       Location_in = T_vals[1, "Location_in"],
                       Net_in = T_vals[1, "Net_in"],
                       Country_in = T_vals[1, "Country_in"],
                       
                       int_Net_in = T_vals[1, "Net_in"],
                       
                       top_up_name = "top_up_IG2_T_0.25y",
                       
                       top_up_net = "IG2",
                       
                       sites_in = T_site,
                       
                       intervention_time = int_time_T,
                       baseline_pred_days = baseline_time_T)


out_p_1 <- sim_forward(start_EIR = T_vals[1, "start_EIR"],
                          
                          Location_in = T_vals[1, "Location_in"],
                          Net_in = T_vals[1, "Net_in"],
                          Country_in = T_vals[1, "Country_in"],
                          
                          int_Net_in = T_vals[1, "Net_in"],
                          
                          top_up_name = "top_up_IG2_T_1y",
                          
                          top_up_net = "Pyrethroid_only",
                          
                          sites_in = T_site,
                          
                          intervention_time = int_time_T,
                          baseline_pred_days = baseline_time_T)




out_IG2_1 <- sim_forward(start_EIR = T_vals[1, "start_EIR"],
                            
                            Location_in = T_vals[1, "Location_in"],
                            Net_in = T_vals[1, "Net_in"],
                            Country_in = T_vals[1, "Country_in"],
                            
                            int_Net_in = T_vals[1, "Net_in"],
                            
                            top_up_name = "top_up_IG2_T_1y",
                            
                            top_up_net = "IG2",
                            
                            sites_in = T_site,
                            
                            intervention_time = int_time_T,
                            baseline_pred_days = baseline_time_T)

out_p_3 <- sim_forward(start_EIR = T_vals[1, "start_EIR"],
                       
                       Location_in = T_vals[1, "Location_in"],
                       Net_in = T_vals[1, "Net_in"],
                       Country_in = T_vals[1, "Country_in"],
                       
                       int_Net_in = T_vals[1, "Net_in"],
                       
                       top_up_name = "top_up_IG2_T_3y",
                       
                       top_up_net = "Pyrethroid_only",
                       
                       sites_in = T_site,
                       
                       intervention_time = int_time_T,
                       baseline_pred_days = baseline_time_T)




out_IG2_3 <- sim_forward(start_EIR = T_vals[1, "start_EIR"],
                         
                         Location_in = T_vals[1, "Location_in"],
                         Net_in = T_vals[1, "Net_in"],
                         Country_in = T_vals[1, "Country_in"],
                         
                         int_Net_in = T_vals[1, "Net_in"],
                         
                         top_up_name = "top_up_IG2_T_3y",
                         
                         top_up_net = "IG2",
                         
                         sites_in = T_site,
                         
                         intervention_time = int_time_T,
                         baseline_pred_days = baseline_time_T)


plot_grid(
  ggplot() +
    geom_line(data = data.frame(t = out_p_0.25$timestep,
                                prev = out_p_0.25$n_detect_182.5_5110/out_p_0.25$n_182.5_5110),
              aes(x = t, y = prev), col = "skyblue") +
    # geom_line(data = data.frame(t = out_none$timestep,
    #                             prev = out_none$n_detect_182.5_5110/out_none$n_182.5_5110),
    #           aes(x = t, y = prev), col = "darkgreen") + 
    geom_vline(xintercept = int_time_T, linetype = 2) +
    geom_line(data = data.frame(t = out_IG2_0.25$timestep,
                                prev = out_IG2_0.25$n_detect_182.5_5110/out_IG2_0.25$n_182.5_5110),
              aes(x = t, y = prev), col = "orange") +
    coord_cartesian(xlim = c(baseline_time_B - 365, baseline_time_B + 365*50)) +
    theme_bw() +
    ggtitle("Quarterly replacement of nets"),
  
  ggplot() +
  geom_line(data = data.frame(t = out_p_1$timestep,
                              prev = out_p_1$n_detect_182.5_5110/out_p_1$n_182.5_5110),
            aes(x = t, y = prev), col = "skyblue") +
  # geom_line(data = data.frame(t = out_none$timestep,
  #                             prev = out_none$n_detect_182.5_5110/out_none$n_182.5_5110),
  #           aes(x = t, y = prev), col = "darkgreen") + 
  geom_vline(xintercept = int_time_T, linetype = 2) +
  geom_line(data = data.frame(t = out_IG2_1$timestep,
                              prev = out_IG2_1$n_detect_182.5_5110/out_IG2_1$n_182.5_5110),
            aes(x = t, y = prev), col = "orange") +
  coord_cartesian(xlim = c(baseline_time_B - 365, baseline_time_B + 365*50)) +
  theme_bw() +
    ggtitle("Yearly replacement of nets"),
  
  ggplot() +
    geom_line(data = data.frame(t = out_p_3$timestep,
                                prev = out_p_3$n_detect_182.5_5110/out_p_3$n_182.5_5110),
              aes(x = t, y = prev), col = "skyblue") +
    # geom_line(data = data.frame(t = out_none$timestep,
    #                             prev = out_none$n_detect_182.5_5110/out_none$n_182.5_5110),
    #           aes(x = t, y = prev), col = "darkgreen") + 
    geom_vline(xintercept = int_time_T, linetype = 2) +
    geom_line(data = data.frame(t = out_IG2_3$timestep,
                                prev = out_IG2_3$n_detect_182.5_5110/out_IG2_3$n_182.5_5110),
              aes(x = t, y = prev), col = "orange") +
    coord_cartesian(xlim = c(baseline_time_B - 365, baseline_time_B + 365*50)) +
    theme_bw() +
    ggtitle("Replacement of nets every three years")
  
  
)




ggplot() +
  geom_line(data = data.frame(t = out_p_0.25$timestep,
                              cov = out_p_0.25$n_use_net),
            aes(x = t, y = cov), col = "skyblue") +
  # geom_line(data = data.frame(t = out_none$timestep,
  #                             prev = out_none$n_detect_182.5_5110/out_none$n_182.5_5110),
  #           aes(x = t, y = prev), col = "darkgreen") + 
  geom_vline(xintercept = int_time_T, linetype = 2) +
  geom_line(data = data.frame(t = out_IG2_0.25$timestep,
                              cov = out_IG2_0.25$n_use_net),
            aes(x = t, y = cov), col = "orange") +
  coord_cartesian(xlim = c(baseline_time_B - 365, baseline_time_B + 365*10)) +
  theme_bw()

ggplot() +
  geom_line(data = data.frame(t = out_p_1$timestep,
                              EIR = out_p_1$EIR_funestus),
            aes(x = t, y = EIR), col = "skyblue") +
  # geom_line(data = data.frame(t = out_none$timestep,
  #                             prev = out_none$n_detect_182.5_5110/out_none$n_182.5_5110),
  #           aes(x = t, y = prev), col = "darkgreen") + 
  geom_vline(xintercept = int_time_T, linetype = 2) +
  geom_line(data = data.frame(t = out_IG2_1$timestep,
                              EIR = out_IG2_1$EIR_funestus),
            aes(x = t, y = EIR), col = "orange") +
  coord_cartesian(xlim = c(baseline_time_B - 365, baseline_time_B + 365*5)) +
  theme_bw()

ggplot() +
  geom_line(data = data.frame(t = out_p_1$timestep,
                              prev = out_p_1$n_detect_182.5_5110/out_p_1$n_182.5_5110),
            aes(x = t, y = prev), col = "skyblue") +
  # geom_line(data = data.frame(t = out_none$timestep,
  #                             prev = out_none$n_detect_182.5_5110/out_none$n_182.5_5110),
  #           aes(x = t, y = prev), col = "darkgreen") + 
  geom_vline(xintercept = int_time_T, linetype = 2) +
  geom_line(data = data.frame(t = out_IG2_1$timestep,
                              prev = out_IG2_1$n_detect_182.5_5110/out_IG2_1$n_182.5_5110),
            aes(x = t, y = prev), col = "orange") +
  coord_cartesian(xlim = c(baseline_time_B - 365, baseline_time_B + 365*5)) +
  theme_bw()


# running the simulations
T_vals_pred_u <- rbind(T_vals %>% mutate(int_Net_in = "IG2", start_EIR = unlist(start_EIR_T)),
                     T_vals %>% mutate(int_Net_in = "PBO", start_EIR = unlist(start_EIR_T)),
                     T_vals %>% mutate(int_Net_in = "Pyrethroid_only", start_EIR = unlist(start_EIR_T)))

T_vals_pred_u <- subset(T_vals_pred_u, Net_in != "RG")

T_vals_pred_b <- rbind(T_vals_pred_u %>% mutate(bioassay_uncertainty = "middle"),
                       T_vals_pred_u %>% mutate(bioassay_uncertainty = "lower"),
                       T_vals_pred_u %>% mutate(bioassay_uncertainty = "upper")) %>% as.data.frame()

T_vals_pred <- bind_rows(
  lapply(seq(1, 5, 1), function(i, T_vals_pred_b){
    T_vals_pred_b %>% mutate(rep = i)}, 
    T_vals_pred_b = T_vals_pred_b)
  )

cl <- makeCluster(5)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("T_site", "int_time_T", "baseline_time_T", "df", "df_species_in", 
                                   "dat_res_pyr", "sim_length", "dat_res_pp", "dat_res_pbo",
                                   "top_up_IG2_T", "top_up_PBO_T", 
                                   "top_up_RG_T", "top_up_p_only_T"))
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
                        bioassay_uncertainty = T_vals_pred[i, "bioassay_uncertainty"])
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
  lapply(seq(1, 5, 1), function(i, T_vals_pred_b){
    T_vals_pred_b %>% mutate(rep = i)}, 
    T_vals_pred_b = T_vals_pred_tu_s)
)

cl <- makeCluster(5)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("T_site", "int_time_T", "baseline_time_T", "df", "df_species_in", 
                                   "dat_res_pyr", "sim_length", "dat_res_pp", "dat_res_pbo",
                                   "top_up_IG2_T", "top_up_PBO_T", 
                                   "top_up_RG_T", "top_up_p_only_T"))
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
                        bioassay_uncertainty = T_vals_pred_tu[i, "bioassay_uncertainty"])
  },
  error = function(cond){
    return(NA)
  })
}
saveRDS(pred_prev_T_tu, 
        file = "data/pred_prev_T_tu.rds")

stopCluster(cl)

# Benin
B_vals_pred_u <- rbind(B_vals %>% mutate(int_Net_in = "IG2", start_EIR = unlist(start_EIR_B)),
                     B_vals %>% mutate(int_Net_in = "PBO", start_EIR = unlist(start_EIR_B)),
                     B_vals %>% mutate(int_Net_in = "Pyrethroid_only", start_EIR = unlist(start_EIR_B)))

B_vals_pred_u <- subset(B_vals_pred_u, Net_in != "RG")

B_vals_pred_b <- rbind(B_vals_pred_u %>% mutate(bioassay_uncertainty = "middle"),
                       B_vals_pred_u %>% mutate(bioassay_uncertainty = "lower"),
                       B_vals_pred_u %>% mutate(bioassay_uncertainty = "upper")) 

B_vals_pred <- bind_rows(
  lapply(seq(1, 5, 1), function(i, B_vals_pred_b){
    B_vals_pred_b %>% mutate(rep = i)}, 
    B_vals_pred_b = B_vals_pred_b)
)

cl <- makeCluster(5)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("B_site", "int_time_B", "baseline_time_B", "df", "df_species_in", 
                                   "dat_res_pyr", "sim_length", "dat_res_pp", "dat_res_pbo",
                                   "top_up_IG2_B",
                                   "top_up_RG_B", "top_up_p_only_B"))
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
                        bioassay_uncertainty = B_vals_pred[i, "bioassay_uncertainty"])
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
  lapply(seq(1, 5, 1), function(i, T_vals_pred_b){
    T_vals_pred_b %>% mutate(rep = i)}, 
    T_vals_pred_b = B_vals_pred_tu_s)
)

cl <- makeCluster(5)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("B_site", "int_time_B", "baseline_time_B", "df", "df_species_in", 
                                   "dat_res_pyr", "sim_length", "dat_res_pp", "dat_res_pbo",
                                   "top_up_IG2_B",
                                   "top_up_RG_B", "top_up_p_only_B"))
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
                        bioassay_uncertainty = B_vals_pred_tu[i, "bioassay_uncertainty"])
  },
  error = function(cond){
    return(NA)
  })
}
saveRDS(pred_prev_B_tu, 
        file = "data/pred_prev_B_tu.rds")

stopCluster(cl)

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

for(i in 1:length(pred_prev_B_RG)){
  pred_prev_B_RG[[i]]$date <- seq(start_date, by = "day", length.out = length(pred_prev_B_RG[[i]]$timestep))
}

for(i in 1:length(pred_prev_T)){
  pred_prev_T[[i]]$date <- seq(start_date, by = "day", length.out = length(pred_prev_T[[1]]$timestep))
}

extract_mean_prev <- function(Net_in, 
                              int_Net_in,
                              Country_in,
                              bioassay_uncertainty,
                              vals_pred, 
                              pred_prev,
                              top_up_net = NA,
                              n_rep = 5){
  
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
    prev_mat <- if(Country_in == "Tanzania"){sapply(inds, function(i, p){summary_pfpr_0.5_5(p[[i]])}, p = pred_prev)
      } else if(Country_in == "Benin"){sapply(inds, function(i, p){summary_pfpr_0_100_all(p[[i]])}, p = pred_prev)}
    return(rowMeans(prev_mat))
  }
}

arm_sims_T <- data.frame("t" = rep(pred_prev_T[[1]]$timestep, 3),
                         "date" = rep(pred_prev_T[[1]]$date, 3),
                         "net" = c(rep("IG2", sim_length),
                                   rep("PBO", sim_length),
                                   rep("Pyrethroid_only", sim_length)),
                         "prev" = c(extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
                                                      Country_in = "Tanzania", bioassay_uncertainty = "middle",
                                                      vals_pred = T_vals_pred, pred_prev = pred_prev_T),
                                    extract_mean_prev(Net_in = "PBO", int_Net_in = "PBO", 
                                                      Country_in = "Tanzania", bioassay_uncertainty = "middle",
                                                      vals_pred = T_vals_pred, pred_prev = pred_prev_T),
                                    extract_mean_prev(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only", 
                                                      Country_in = "Tanzania", bioassay_uncertainty = "middle",
                                                      vals_pred = T_vals_pred, pred_prev = pred_prev_T)
                                    ),
                         "lower" = c(extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
                                                      Country_in = "Tanzania", bioassay_uncertainty = "lower",
                                                      vals_pred = T_vals_pred, pred_prev = pred_prev_T),
                                    extract_mean_prev(Net_in = "PBO", int_Net_in = "PBO", 
                                                      Country_in = "Tanzania", bioassay_uncertainty = "lower",
                                                      vals_pred = T_vals_pred, pred_prev = pred_prev_T),
                                    extract_mean_prev(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only", 
                                                      Country_in = "Tanzania", bioassay_uncertainty = "lower",
                                                      vals_pred = T_vals_pred, pred_prev = pred_prev_T)
                         ),
                         "upper" = c(extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
                                                      Country_in = "Tanzania", bioassay_uncertainty = "upper",
                                                      vals_pred = T_vals_pred, pred_prev = pred_prev_T),
                                    extract_mean_prev(Net_in = "PBO", int_Net_in = "PBO", 
                                                      Country_in = "Tanzania", bioassay_uncertainty = "upper",
                                                      vals_pred = T_vals_pred, pred_prev = pred_prev_T),
                                    extract_mean_prev(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only", 
                                                      Country_in = "Tanzania", bioassay_uncertainty = "upper",
                                                      vals_pred = T_vals_pred, pred_prev = pred_prev_T)
                         )
                         )


arm_sims_B <- data.frame("t" = rep(pred_prev_B[[1]]$timestep, 2),
                         "date" = rep(pred_prev_B[[1]]$date, 2),
                         "net" = c(rep("IG2", sim_length),
                                   rep("Pyrethroid_only", sim_length)),
                         "prev" = c(extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
                                                      Country_in = "Benin", bioassay_uncertainty = "middle",
                                                      vals_pred = B_vals_pred, pred_prev = pred_prev_B),
                                    extract_mean_prev(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only", 
                                                      Country_in = "Benin", bioassay_uncertainty = "middle",
                                                      vals_pred = B_vals_pred, pred_prev = pred_prev_B)
                         ),
                         "lower" = c(extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
                                                       Country_in = "Benin", bioassay_uncertainty = "lower",
                                                       vals_pred = B_vals_pred, pred_prev = pred_prev_B),
                                     extract_mean_prev(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only", 
                                                       Country_in = "Benin", bioassay_uncertainty = "lower",
                                                       vals_pred = B_vals_pred, pred_prev = pred_prev_B)
                         ),
                         "upper" = c(extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
                                                       Country_in = "Benin", bioassay_uncertainty = "upper",
                                                       vals_pred = B_vals_pred, pred_prev = pred_prev_B),
                                     extract_mean_prev(Net_in = "Pyrethroid_only", int_Net_in = "Pyrethroid_only", 
                                                       Country_in = "Benin", bioassay_uncertainty = "upper",
                                                       vals_pred = B_vals_pred, pred_prev = pred_prev_B)
                         )
)

# Tanzania
T_arm_plot <- ggplot(data = arm_sims_T) +
  geom_vline(xintercept = intervention_date_T, linetype = 2, linewidth = 1, alpha = 0.5) +
  geom_ribbon(aes(x = date, ymin = lower, ymax = upper, fill = net), alpha = 0.25) +
  geom_line(aes(x = date, y = prev, col = net)) +
  geom_pointrange(data = subset(df_T, Net == "PBO"), 
             aes(x = date, y = Malaria_prevalence, ymin = Malaria_prevalence_l, ymax = Malaria_prevalence_u),
             fill = "purple", size = 0.75, shape = 21) +
  geom_pointrange(data = subset(df_T, Net == "IG2"), 
             aes(x = date, y = Malaria_prevalence, ymin = Malaria_prevalence_l, ymax = Malaria_prevalence_u),
             fill = "aquamarine4", size = 0.75, shape = 21) +
  geom_pointrange(data = subset(df_T, Net == "Pyrethroid_only"), 
             aes(x = date, y = Malaria_prevalence, ymin = Malaria_prevalence_l, ymax = Malaria_prevalence_u),
             fill = "black", size = 0.75, shape = 21) +
  scale_colour_manual(values = c("aquamarine4", "purple", "black"), name = "",
                      labels = c("pyrethroid-pyrrole", "pyrethroid-PBO", "pyrethroid-only")) +
  scale_fill_manual(values = c( "aquamarine4","purple", "black"), name = "",
                    labels = c("pyrethroid-pyrrole", "pyrethroid-PBO", "pyrethroid-only")) +
  xlab("Year") +
  ylab("Malaria prevalence in children\naged 0.5 to 14 years old") +
  coord_cartesian(xlim = c(baseline_start_date_B - 365, baseline_start_date_B + 365*5)) +
  theme_classic() +
  ggtitle("Tanzania") +
  theme(text = element_text(size = 17)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), breaks = seq(0, 1, 0.2))

B_arm_plot <- ggplot(data = arm_sims_B) +
  geom_vline(xintercept = intervention_date_B, linetype = 2, linewidth = 1, alpha = 0.5) +
  geom_ribbon(aes(x = date, ymin = lower, ymax = upper, fill = net), alpha = 0.25) +
  geom_line(aes(x = date, y = prev, col = net)) +
  # geom_point(data = subset(df_B, Location == "Cove"), 
  #            aes(x = date, y = Malaria_prevalence),
  #            fill = "skyblue", size = 3.5, shape = 21) +
  geom_pointrange(data = subset(df_B, Location == "Zagnanado"), 
                  aes(x = date, y = Malaria_prevalence, ymin = Malaria_prevalence_l, ymax = Malaria_prevalence_u),
                  fill = "aquamarine4", size = 0.75, shape = 22) +
  geom_pointrange(data = subset(df_B, Location == "Ouinhi"), 
                  aes(x = date, y = Malaria_prevalence, ymin = Malaria_prevalence_l, ymax = Malaria_prevalence_u),
                  fill = "black", size = 0.75, shape = 22) +
  scale_colour_manual(values = c("aquamarine4", "black"), name = "") +
  scale_fill_manual(values = c("aquamarine4", "black"), name = "") +
  xlab("Year") +
  ylab("Malaria prevalence in\npeople of all ages") +
  coord_cartesian(xlim = c(baseline_start_date_B - 365, baseline_start_date_B + 365*3)) +
  theme_classic() +
  ggtitle("Benin") +
  theme(text = element_text(size = 17)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), breaks = seq(0, 1, 0.2))

# pyrethroid-pyrrole arms with pyrethroid-pyrrole top up nets
T_tu_df <- data.frame("t" = rep(pred_prev_T[[1]]$timestep, 3),
                      "date" = rep(pred_prev_T[[1]]$date, 3),
                      "net" = c(rep("IG2", sim_length),
                                rep("IG2_tu", sim_length),
                                rep("Pyrethroid_cf", sim_length)),
                      "prev" = c(extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
                                                   Country_in = "Tanzania", bioassay_uncertainty = "middle",
                                                   vals_pred = T_vals_pred, pred_prev = pred_prev_T),
                                 
                                 extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
                                                   Country_in = "Tanzania", bioassay_uncertainty = "middle",
                                                   vals_pred = T_vals_pred_tu, pred_prev = pred_prev_T_tu,
                                                   top_up_net = "IG2"),
                                 
                                 extract_mean_prev(Net_in = "IG2", int_Net_in = "Pyrethroid_only", 
                                                   Country_in = "Tanzania", bioassay_uncertainty = "middle",
                                                   vals_pred = T_vals_pred, pred_prev = pred_prev_T)
                      ),
                      "lower" = c(extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
                                                    Country_in = "Tanzania", bioassay_uncertainty = "lower",
                                                    vals_pred = T_vals_pred, pred_prev = pred_prev_T),
                                  
                                  extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
                                                    Country_in = "Tanzania", bioassay_uncertainty = "lower",
                                                    vals_pred = T_vals_pred_tu, pred_prev = pred_prev_T_tu,
                                                    top_up_net = "IG2"),
                                  
                                  extract_mean_prev(Net_in = "IG2", int_Net_in = "Pyrethroid_only", 
                                                    Country_in = "Tanzania", bioassay_uncertainty = "lower",
                                                    vals_pred = T_vals_pred, pred_prev = pred_prev_T)
                      ),
                      "upper" = c(extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
                                                    Country_in = "Tanzania", bioassay_uncertainty = "upper",
                                                    vals_pred = T_vals_pred, pred_prev = pred_prev_T),
                                  
                                  extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
                                                    Country_in = "Tanzania", bioassay_uncertainty = "upper",
                                                    vals_pred = T_vals_pred_tu, pred_prev = pred_prev_T_tu,
                                                    top_up_net = "IG2"),
                                  
                                  extract_mean_prev(Net_in = "IG2", int_Net_in = "Pyrethroid_only", 
                                                    Country_in = "Tanzania", bioassay_uncertainty = "upper",
                                                    vals_pred = T_vals_pred, pred_prev = pred_prev_T)
                      )
)

ggplot(data = T_tu_df, 
       aes(x = date, y = prev, ymin = lower, ymax = upper, fill = net)) +
  geom_ribbon(alpha = 0.25) +
  geom_line(aes(col = net)) +
  coord_cartesian(xlim = c(baseline_start_date_B - 365, baseline_start_date_B + 365*5))

# pyrethroid-pyrrole arms with pyrethroid-pyrrole top up nets
B_tu_df <- data.frame("t" = rep(pred_prev_B[[1]]$timestep, 3),
                      "date" = rep(pred_prev_B[[1]]$date, 3),
                      "net" = c(rep("IG2", sim_length),
                                rep("IG2_tu", sim_length),
                                rep("Pyrethroid_cf", sim_length)),
                      "prev" = c(extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
                                                   Country_in = "Benin", bioassay_uncertainty = "middle",
                                                   vals_pred = B_vals_pred, pred_prev = pred_prev_B),
                                 
                                 extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
                                                   Country_in = "Benin", bioassay_uncertainty = "middle",
                                                   vals_pred = B_vals_pred_tu, pred_prev = pred_prev_B_tu,
                                                   top_up_net = "IG2"),
                                 
                                 extract_mean_prev(Net_in = "IG2", int_Net_in = "Pyrethroid_only", 
                                                   Country_in = "Benin", bioassay_uncertainty = "middle",
                                                   vals_pred = B_vals_pred, pred_prev = pred_prev_B)
                      ),
                      "lower" = c(extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
                                                    Country_in = "Benin", bioassay_uncertainty = "lower",
                                                    vals_pred = B_vals_pred, pred_prev = pred_prev_B),
                                  
                                  extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
                                                    Country_in = "Benin", bioassay_uncertainty = "lower",
                                                    vals_pred = B_vals_pred_tu, pred_prev = pred_prev_B_tu,
                                                    top_up_net = "IG2"),
                                  
                                  extract_mean_prev(Net_in = "IG2", int_Net_in = "Pyrethroid_only", 
                                                    Country_in = "Benin", bioassay_uncertainty = "lower",
                                                    vals_pred = B_vals_pred, pred_prev = pred_prev_B)
                      ),
                      "upper" = c(extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
                                                    Country_in = "Benin", bioassay_uncertainty = "upper",
                                                    vals_pred = B_vals_pred, pred_prev = pred_prev_B),
                                  
                                  extract_mean_prev(Net_in = "IG2", int_Net_in = "IG2", 
                                                    Country_in = "Benin", bioassay_uncertainty = "upper",
                                                    vals_pred = B_vals_pred_tu, pred_prev = pred_prev_B_tu,
                                                    top_up_net = "IG2"),
                                  
                                  extract_mean_prev(Net_in = "IG2", int_Net_in = "Pyrethroid_only", 
                                                    Country_in = "Benin", bioassay_uncertainty = "upper",
                                                    vals_pred = B_vals_pred, pred_prev = pred_prev_B)
                      )
)

ggplot(data = B_tu_df, 
       aes(x = date, y = prev, ymin = lower, ymax = upper, fill = net)) +
  geom_ribbon(alpha = 0.25) +
  geom_line(aes(col = net)) +
  coord_cartesian(xlim = c(baseline_start_date_B - 365, baseline_start_date_B + 365*3))

# actual vs fitted plot
inds_T <- match(interaction(df_T$date, df_T$Net), interaction(arm_sims_T$date, arm_sims_T$net))

df_T <- df_T %>% mutate(pred_prev = arm_sims_T[inds_T, "prev"],
                        pred_prev_l = arm_sims_T[inds_T, "lower"],
                        pred_prev_u = arm_sims_T[inds_T, "upper"])

inds_B <- match(interaction(df_B$date, df_B$Net), interaction(arm_sims_B$date, arm_sims_B$net))

df_B <- df_B %>% mutate(pred_prev = arm_sims_B[inds_B, "prev"],
                        pred_prev_l = arm_sims_B[inds_B, "lower"],
                        pred_prev_u = arm_sims_B[inds_B, "upper"])

df_af <- rbind(df_T[, c("Net", "Country", "date", "Malaria_prevalence", "Malaria_prevalence_u", 
                        "Malaria_prevalence_l", "pred_prev", "pred_prev_l", "pred_prev_u")], 
               df_B[, c("Net", "Country", "date", "Malaria_prevalence", "Malaria_prevalence_u",
                        "Malaria_prevalence_l", "pred_prev", "pred_prev_l", "pred_prev_u")])

df_af <- df_af[is.na(df_af$Malaria_prevalence)!=1, ]

af_plot <- ggplot(data = subset(df_af, !(date %in% c(baseline_start_date_B,
                                                     baseline_start_date_T)) & Net != "RG")) +
  geom_abline(slope = 1, linetype = 2, linewidth = 0.75, alpha = 0.4) +
  #geom_smooth(formula = y ~ x, se = FALSE, aes(x = Malaria_prevalence, y = pred_prev),
  #            method = "lm", col = "grey50", fullrange = TRUE, alpha = 0.8) +
  geom_pointrange(size = 0.75, 
                  aes(x = Malaria_prevalence, y = pred_prev, ymin = pred_prev_l, ymax = pred_prev_u, col = Net, shape = Country)) +
  geom_errorbarh(aes(y = pred_prev, xmin = Malaria_prevalence_l, xmax = Malaria_prevalence_u, col = Net, group = Country)) +
  theme_classic() + theme(text = element_text(size = 18)) +
  scale_colour_manual(values = c("aquamarine4", "purple", "black"),
                      labels = c("pyrethroid-pyrrole", "pyrethroid-PBO", "pyrethroid-only")) +
  scale_shape_manual(values = c(15, 16)) +
  ylab("Predicted malaria prevalence") + xlab("Observed malaria prevalence") +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent, breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent, breaks = seq(0, 1, 0.2))

plot_grid(
plot_grid(
  T_arm_plot + theme(legend.position = "none"),
  
  B_arm_plot + 
    theme(legend.position = "none"),
  
  #get_legend(T_arm_plot),
  labels = c("A", "B", ""), nrow = 1#, rel_widths = c(1, 1, 0.4)
),
plot_grid(af_plot, 
  NULL, labels = c("C", ""), rel_widths = c(1, 0.5),
  nrow = 1), nrow = 2
)


  
# calculate efficacy - put on a different plot
# put grey points from previous paper - R2 with and without the new ones
# R2 from the x=y line

# Joe has measurement error

# Counterfactuals
# IG2

calc_incidence <- function(output, s_time, e_time){
  sum(output$n_inc_clinical_0_36500[s_time:e_time]/output$n_0_36500[s_time:e_time])
}

p_only_1y <- calc_incidence(out_pyrethroid_only, intervention_time, intervention_time + 365)
p_only_2y <- calc_incidence(out_pyrethroid_only, intervention_time, intervention_time + 365*2)
ppf_1y <- calc_incidence(out_ppf, intervention_time, intervention_time + 365)
ppf_2y <- calc_incidence(out_ppf, intervention_time, intervention_time + 365*2)
pp_1y <- calc_incidence(out_IG2, intervention_time, intervention_time + 365)
pp_2y <- calc_incidence(out_IG2, intervention_time, intervention_time + 365*2)

rel_reduction <- rbind(data.frame(year = c("1", "2"),
                                  ITN = c("pyrrole", "pyrrole"),
                                  r = c(((p_only_1y - pp_1y)/p_only_1y)*100,
                                        ((p_only_2y - pp_2y)/p_only_2y)*100
                                  )),
                       data.frame(year = c("1", "2"),
                                  ITN = c("pyriproxyfen", "pyriproxyfen"),
                                  r = c(((p_only_1y - ppf_1y)/p_only_1y)*100,
                                        ((p_only_2y - ppf_2y)/p_only_2y)*100
                                  ))
)

ggplot(data = rel_reduction,
       aes(x = year, y = r, col = ITN)) +
  geom_point(size = 3) + theme_classic() +
  ylab("Relative reduction in incidence") + xlab("Year") +
  scale_colour_manual(values = c("skyblue", "aquamarine4"),
                      labels = c("Pyrethroid-pyriproxyfen", "Pyrethroid-pyrrole")) +
  ggtitle("Benin") +
  scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100)) +
  theme(text = element_text(size = 18))

write.csv(data.frame("t" = out_IG2$timestep,
                     "date" = out_IG2$date,
                     "RG" = summary_pfpr_0_100_all(pred_prev_B_RG[[2]]),
                     "IG2" = summary_pfpr_0_100_all(pred_prev_B[[2]]),
                     "pbo" = summary_pfpr_0_100_all(pred_prev_B[[5]]),
                     "pyr" = summary_pfpr_0_100_all(pred_prev_B[[8]])),
          file = "IG2_arm_cf.csv")

ggplot(data = data.frame("t" = out_IG2$timestep,
                         "date" = out_IG2$date,
                         "RG" = summary_pfpr_0_100_all(pred_prev_B_RG[[2]]),
                         "IG2" = summary_pfpr_0_100_all(pred_prev_B[[2]]),
                         "pbo" = summary_pfpr_0_100_all(pred_prev_B[[5]]),
                         "pyr" = summary_pfpr_0_100_all(pred_prev_B[[8]]))) +
  geom_line(aes(x = date, y = IG2), col = "aquamarine4") +
  geom_line(aes(x = date, y = RG), col = "skyblue") +
  geom_line(aes(x = date, y = pbo), col = "purple") +
  geom_line(aes(x = date, y = pyr)) + 
  geom_point(data = subset(df_B, Location == "Zagnanado"), 
             aes(x = date, y = Malaria_prevalence),
             fill = "aquamarine4", size = 3.5, shape = 21) +
  xlab("Year") +
  ylab("Malaria prevalence in total population") +
  coord_cartesian(xlim = c(baseline_start_date_B - 365, baseline_start_date_B + 365*3)) +
  theme_classic() + ylim(0, 1)

#pyrethroid_only
write.csv(data.frame("t" = out_IG2$timestep,
                     "date" = out_IG2$date,
                     "RG" = summary_pfpr_0_100_all(pred_prev_B_RG[[3]]),
                     "IG2" = summary_pfpr_0_100_all(pred_prev_B[[3]]),
                     "pbo" = summary_pfpr_0_100_all(pred_prev_B[[6]]),
                     "pyr" = summary_pfpr_0_100_all(pred_prev_B[[9]])),
          file = "p_only_arm_cf.csv")

ggplot(data = data.frame("t" = out_IG2$timestep,
                         "date" = out_IG2$date,
                         "RG" = summary_pfpr_0_100_all(pred_prev_B_RG[[3]]),
                         "IG2" = summary_pfpr_0_100_all(pred_prev_B[[3]]),
                         "pbo" = summary_pfpr_0_100_all(pred_prev_B[[6]]),
                         "pyr" = summary_pfpr_0_100_all(pred_prev_B[[9]]))) +
  geom_line(aes(x = date, y = IG2), col = "aquamarine4") +
  geom_line(aes(x = date, y = RG), col = "skyblue") +
  geom_line(aes(x = date, y = pbo), col = "purple") +
  geom_line(aes(x = date, y = pyr)) +
  geom_point(data = subset(df_B, Location == "Ouinhi"), 
             aes(x = date, y = Malaria_prevalence),
             fill = "black", size = 3.5, shape = 21) +
  ylab("Malaria prevalence in total population") +
  coord_cartesian(xlim = c(baseline_start_date_B - 365, baseline_start_date_B + 365*3)) +
  theme_classic() + ylim(0, 1)

# RG

write.csv(data.frame("t" = out_IG2$timestep,
                     "date" = out_IG2$date,
                     "RG" = summary_pfpr_0_100_all(pred_prev_B_RG[[1]]),
                     "IG2" = summary_pfpr_0_100_all(pred_prev_B[[1]]),
                     "pbo" = summary_pfpr_0_100_all(pred_prev_B[[4]]),
                     "pyr" = summary_pfpr_0_100_all(pred_prev_B[[7]])),
          file = "RG_arm_cf.csv")

saveRDS(out_IG2, file = "out_IG2_B.rds")
saveRDS(out_ppf, file = "out_ppf_B.rds")
saveRDS(out_pyrethroid_only, file = "out_pyrethroid_only_B.rds")

pred_times <- unique(c(df$Time_months, 36)*30) + int_time
pred_times[1] <- pred_times[1] - 1

pred_indices <- which(out_IG2$timestep %in% pred_times)

round(summary_pfpr_0_100(out_IG2)[pred_indices], digits = 2)

round(summary_pfpr_0_100(out_pyrethroid_only)[pred_indices], digits = 2)

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

df <- subset(df, Country == "Benin")

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