suppressPackageStartupMessages(library(ggplot2)); library(malariasimulation); library(malariaEquilibrium)
library(reshape2); library(tidyverse); library(readxl); library(rstan); library(pracma)
library(foresite); library(doParallel); library(foreach);
library(patchwork); library(pammtools); library(ggpattern);
library(cowplot)

####################################
##### Read in seasonality data #####
####################################
# getting the site files for Benin and Tanzania

B_site <- foresite::BEN
T_site <- foresite::TZA

##########################################################
##### RCT nets retention and malaria prevalence data #####
##########################################################

# Tanzania: https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#supplementaryMaterial
# Benin: https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(22)02319-4/fulltext#supplementaryMaterial

df <- read_excel("data/data_RCT.xlsx", sheet = "df") %>% 
  mutate(Bed_net_use_both = Bed_net_use_both/100,
         Bed_net_use_RCT = Bed_net_use_RCT/100,
         Malaria_prevalence = Malaria_prevalence/100,
         Malaria_prevalence_l = Malaria_prevalence - 1.96*sqrt(Malaria_prevalence * (1 - Malaria_prevalence) / Malaria_prev_n_tested),
         Malaria_prevalence_u = Malaria_prevalence + 1.96*sqrt(Malaria_prevalence * (1 - Malaria_prevalence) / Malaria_prev_n_tested))

df_species <- read_excel("data/data_RCT.xlsx", sheet = "site_params")

#############################################################
##### read in the parameter efficacy estimates for nets #####
#############################################################

dat_res_pyr_ll <- readRDS("parameters/pyrethroid_uncertainty 2") %>% dplyr::rename(rn0 = rn_pyr,
                                                                           gamman = mean_duration) %>% 
  dplyr::mutate(bioassay_mortality = 1 - resistance)  %>% select(dn0, rn0, gamman, resistance, bioassay_mortality)

dat_res_pbo_ll <- readRDS("parameters/pbo_uncertainty_using_pyrethroid_dn0_for_mn_durability 2") %>% dplyr::rename(rn0 = rn_pbo,
                                                                                                           gamman = mn_dur,
                                                                                                           resistance = resistance.x) %>%
  
  dplyr::mutate(bioassay_mortality = 1 - resistance) %>% select(dn0, rn0, gamman, resistance, bioassay_mortality)

dat_res_pp_ll <- readRDS("parameters/pyrrole_uncertainty_using_pyrethroid_dn0_for_mn_durability 2") %>% dplyr::rename(rn0 = rn_pbo,
                                                                                                              gamman = mn_dur,
                                                                                                              resistance = resistance.x) %>%
  
  dplyr::mutate(bioassay_mortality = 1 - resistance) %>% select(dn0, rn0, gamman, resistance, bioassay_mortality)


dat_res_pyr <- readRDS("parameters/pyrethroid_binomial_uncertainty 1") %>% dplyr::rename(rn0 = rn_pyr,
                                                                                     gamman = mean_duration) %>% 
  dplyr::mutate(bioassay_mortality = 1 - resistance)  %>% select(dn0, rn0, gamman, resistance, bioassay_mortality)

dat_res_pbo <- readRDS("parameters/pbo_binomial_uncertainty_using_pyrethroid_dn0_for_mn_durability 1") %>% dplyr::rename(rn0 = rn_pbo,
                                                                                                                     gamman = mn_dur,
                                                                                                                     resistance = resistance.x) %>%
  
  dplyr::mutate(bioassay_mortality = 1 - resistance) %>% select(dn0, rn0, gamman, resistance, bioassay_mortality)

dat_res_pp <- readRDS("parameters/pyrrole_binomial_uncertainty_using_pyrethroid_dn0_for_mn_durability 1") %>% dplyr::rename(rn0 = rn_pbo,
                                                                                                                        gamman = mn_dur,
                                                                                                                        resistance = resistance.x) %>%
  
  dplyr::mutate(bioassay_mortality = 1 - resistance) %>% select(dn0, rn0, gamman, resistance, bioassay_mortality)

# functions to generate all the required values
gen_inds <- function(df_w){
  return(df_w %>% 
           dplyr::group_by(resistance, bioassay_mortality) %>% 
           dplyr::mutate(index = row_number()) %>%
           ungroup() %>%
           as.data.frame()
         )
}

gen_median <- function(df_w){
  return(df_w %>% dplyr::group_by(resistance, bioassay_mortality) %>% 
           dplyr::summarise(dn0_med = median(dn0),
                            rn0_med = median(rn0),
                            gamman_med = median(gamman)
                            ) %>%
           ungroup() %>%
           as.data.frame()
  )
}

dat_res_pyr_ll_med <- gen_median(dat_res_pyr_ll)
dat_res_pyr_ll <- gen_inds(dat_res_pyr_ll) 

dat_res_pbo_ll_med <- gen_median(dat_res_pbo_ll)
dat_res_pbo_ll <- gen_inds(dat_res_pbo_ll) 

dat_res_pp_ll_med <- gen_median(dat_res_pp_ll)
dat_res_pp_ll <- gen_inds(dat_res_pp_ll) 

dat_res_pyr_med <- gen_median(dat_res_pyr)
dat_res_pyr <- gen_inds(dat_res_pyr) 

dat_res_pbo_med <- gen_median(dat_res_pbo)
dat_res_pbo <- gen_inds(dat_res_pbo) 

dat_res_pp_med <- gen_median(dat_res_pp)
dat_res_pp <- gen_inds(dat_res_pp) 

# binomial uncertainty estimates
# dat_res_pyr_both <- read.csv("parameters/loglogistic_beta_binomial_and_binomial_pyrethroid_only_nets.csv",header=TRUE) %>% mutate(resistance = 1 - bioassay_mortality)
# dat_res_pbo_both <- read.csv("parameters/loglogistic_beta_binomial_and_binomial_pyrethroid_pbo_nets.csv",header=TRUE) %>% mutate(resistance = 1 - bioassay_mortality)
# dat_res_pp_both <- read.csv("parameters/loglogistic_beta_binomial_and_binomial_pyrethroid_pyrrole_nets.csv",header=TRUE) %>% mutate(resistance = 1 - bioassay_mortality)
# 
# dat_res_pyr <- dat_res_pyr_both %>% rename(dn0_lo = dn0_binomial_05, dn0_up = dn0_binomial_up95, rn0_lo = rn0_binomial_05, rn0_up = rn0_binomial_up95, 
#                                            gamman_lo = gamman_binomial_lo05, gamman_up = gamman_binomial_up95) %>% select(dn0_lo, dn0_med, dn0_up, rn0_lo, rn0_med, rn0_up, gamman_lo, gamman_med, gamman_up,
#                                                                                                                           bioassay_mortality, resistance)
# 
# dat_res_pbo <- dat_res_pbo_both %>% rename(dn0_lo = dn0_binomial_05, dn0_up = dn0_binomial_up95, rn0_lo = rn0_binomial_05, rn0_up = rn0_binomial_up95, 
#                                            gamman_lo = gamman_binomial_lo05, gamman_up = gamman_binomial_up95) %>% select(dn0_lo, dn0_med, dn0_up, rn0_lo, rn0_med, rn0_up, gamman_lo, gamman_med, gamman_up,
#                                                                                                                           bioassay_mortality, resistance)
# dat_res_pp <- dat_res_pp_both %>% rename(dn0_lo = dn0_binomial_05, dn0_up = dn0_binomial_up95, rn0_lo = rn0_binomial_05, rn0_up = rn0_binomial_up95, 
#                                          gamman_lo = gamman_binomial_lo05, gamman_up = gamman_binomial_up95) %>% select(dn0_lo, dn0_med, dn0_up, rn0_lo, rn0_med, rn0_up, gamman_lo, gamman_med, gamman_up,
#                                                                                                                         bioassay_mortality, resistance)
# 
# dat_res_pyr_ll <- dat_res_pyr_both %>% rename(dn0_lo = dn0_betabin_lo05, dn0_up = dn0_betabin_up95, rn0_lo = rn0_betabin_lo05, rn0_up = rn0_betabin_up95, 
#                                               gamman_lo = gamman_betabin_lo05, gamman_up = gamman_betabinom_up95) %>% select(dn0_lo, dn0_med, dn0_up, rn0_lo, rn0_med, rn0_up, gamman_lo, gamman_med, gamman_up,
#                                                                                                                              bioassay_mortality, resistance)
# 
# dat_res_pbo_ll <- dat_res_pbo_both %>% rename(dn0_lo = dn0_betabin_lo05, dn0_up = dn0_betabin_up95, rn0_lo = rn0_betabin_lo05, rn0_up = rn0_betabin_up95, 
#                                               gamman_lo = gamman_betabin_lo05, gamman_up = gamman_betabinom_up95) %>% select(dn0_lo, dn0_med, dn0_up, rn0_lo, rn0_med, rn0_up, gamman_lo, gamman_med, gamman_up,
#                                                                                                                              bioassay_mortality, resistance)
# 
# dat_res_pp_ll <- dat_res_pp_both %>% rename(dn0_lo = dn0_betabin_lo05, dn0_up = dn0_betabin_up95, rn0_lo = rn0_betabin_lo05, rn0_up = rn0_betabin_up95, 
#                                             gamman_lo = gamman_betabin_lo05, gamman_up = gamman_betabinom_up95) %>% select(dn0_lo, dn0_med, dn0_up, rn0_lo, rn0_med, rn0_up, gamman_lo, gamman_med, gamman_up,
#                                                                                                                            bioassay_mortality, resistance)


############################
##### vector bionomics #####
############################

# proportion of bites indoors and in bed from Sherrard-Smith paper (https://www.pnas.org/doi/abs/10.1073/pnas.1820646116)
bite_params <- read_excel("parameters/PNAS.biting.time.xlsx", sheet = "All data for GLMMs")

# calculating the median parameter values by country
bite_params_df <- bite_params %>% group_by(Country_clean, species) %>% summarise(phi_indoor = median(Indoor_phi_mnHuman),
                                                                                 phi_bed = median(Inbed_phi_mnHuman)) %>% 
  mutate(country = ifelse(Country_clean == "BK", "Burkina_Faso", Country_clean))

# top up estimates

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

df$Time_months <- suppressWarnings(as.numeric(df$Time_months)) 

df <- rbind(subset(df, Country == "Tanzania") %>% rowwise() %>% mutate(date = if(is.na(Time_months) == 1){baseline_start_date_T}else{intervention_date_T + months(Time_months)}),
            subset(df, Country == "Benin") %>% rowwise() %>% mutate(date = if(is.na(Time_months) == 1){baseline_start_date_B}else{intervention_date_B + months(Time_months)})
)

df_B <- subset(df, Country == "Benin")

df_T <- subset(df, Country == "Tanzania")

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

