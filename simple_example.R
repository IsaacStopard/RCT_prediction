# simple example of replacing nets with high d0 with low d0
# runs the malaria simulation model where all nets initially implemented have a high or low death probability
# in both scenarios nets are replaced every 10 days with new nets with a low death probability
# all nets are assumed to be retained for a very long time and have a very long half-life

rm(list = ls())

library(malariasimulation); library(malariaEquilibrium); library(tidyverse)

###########################
##### values to input #####
###########################
year <- 365
sim_length <- 50 * year
human_population <- 1000
starting_EIR <- 30

int_time <- 40 * year # time when nets are first added
d0_max <- 0.8 # bed net death probability high value
d0_min <- 0.1 # bed net death probability low value

init_cov <- 0.8 # initial bed net coverage when nets are first implemented

t_diff <- 10
rep_times_trial <- seq(int_time + t_diff, sim_length, t_diff) # times to replace initial nets with (every 10 days)

rep_times_p_only <- rep_times_trial - 1

rep_times <- sort(c(rep_times_p_only, rep_times_trial))

prop_rep <- 0.8 * exp(-0.001 * (rep_times_trial-int_time)) # coverage to replace initial nets with

n_rep_times <- length(rep_times)

covs <- c(rbind(rep(init_cov, length(rep_times_p_only)), prop_rep))
net_names <- c(rbind(paste0("Replacement nets: ", seq(1, length(rep_times_p_only))), paste0("Trial nets: ", seq(1, length(rep_times_trial)))))

d0_params_trial <- c(rbind(rep(d0_min, length(rep_times_p_only)), rep(d0_max, length(rep_times_trial))))

d0_params_no_trial <- rep(d0_min, length(rep_times)) 

###################################
##### running the simulations #####
###################################

simparams <- get_parameters(
  list(
    human_population = human_population,
    model_seasonality = FALSE, # assumes no seasonality
    incidence_rendering_min_ages = 2*year,
    incidence_rendering_max_ages = 10*year
  )
)

# one vector species
simparams <- set_species(simparams,
                         list(gamb_params),
                         c(1))


simparams <- set_equilibrium(simparams, starting_EIR)

bednetparams <- simparams

bednet_events = data.frame(
  timestep = rep_times,
  name= net_names)

# parameters where there are initially more efficacious nets
bednetparams_trial <- set_bednets(
  bednetparams,
  timesteps = bednet_events$timestep,
  coverages = covs,
  retention = 10^10,
  dn0 = matrix(d0_params_trial, nrow = n_rep_times, ncol=1), # replacing the initial trial nets with less efficacious nets
  rn = matrix(rep(0.1, n_rep_times), nrow = n_rep_times, ncol=1), # all other net parameters assumed to remain the same
  rnm = matrix(rep(.24, n_rep_times), nrow = n_rep_times, ncol=1),
  gamman = rep(10^10, n_rep_times) # long bed net half life
)

# parameters where all nets have low efficacy
bednetparams_no_trial <- set_bednets(
  bednetparams,
  timesteps = bednet_events$timestep,
  coverages = rep(init_cov, n_rep_times),
  retention = 10^10,
  dn0 = matrix(d0_params_no_trial, nrow = n_rep_times, ncol=1),
  rn = matrix(rep(0.1, n_rep_times), nrow = n_rep_times, ncol=1),
  rnm = matrix(rep(.24, n_rep_times), nrow = n_rep_times, ncol=1),
  gamman = rep(10^10, n_rep_times)
)

# setting the correlation parameter so bed nets added replace existing bed nets first
correlationsb1 <- get_correlation_parameters(bednetparams_trial)
correlationsb1$inter_round_rho('bednets', 1)

correlationsb2 <- get_correlation_parameters(bednetparams_no_trial)
correlationsb2$inter_round_rho('bednets', 1)

output_trial <- run_simulation(sim_length, bednetparams_trial, correlations = correlationsb1)

output_no_trial <- run_simulation(sim_length, bednetparams_no_trial, correlations = correlationsb2)

#################
##### plots #####
#################

msim_res <- rbind(data.frame(t = output_no_trial$timestep,
                     prev = output_no_trial$n_detect_730_3650/output_no_trial$n_730_3650,
                     EIR = output_no_trial$EIR_gamb,
                     net_use = output_no_trial$n_use_net) %>% mutate(net = "initial nets: low d0"),
                data = data.frame(t = output_trial$timestep,
                              prev = output_trial$n_detect_730_3650/output_trial$n_730_3650,
                              EIR = output_trial$EIR_gamb,
                              net_use = output_trial$n_use_net) %>% mutate(net = "initial nets: high d0"))

# prevalence 
ggplot(data = msim_res) +
  geom_line(aes(x = t, y = prev, col = net))  +
  ylab("Predicted malaria prevalence among 2 to 10 years olds") +
  coord_cartesian(xlim = c(int_time, int_time + 365*10)) + theme_bw()  +
  xlab("Days") +
  scale_colour_manual(values = c("skyblue", "orange"))

# EIR
ggplot(data = msim_res) +
  geom_line(aes(x = t, y = EIR, col = net)) +
  ylab("Simulated EIR") + coord_cartesian(xlim = c(int_time, int_time + 365*70)) + theme_bw()  +
  xlab("Days") +
  scale_colour_manual(values = c("skyblue", "orange"))

# net use
ggplot(data = msim_res) +
  geom_line(aes(x = t, y = net_use, col = net)) +
  ylab("Simulated net use") + coord_cartesian(xlim = c(int_time, int_time + 365*70)) + theme_bw()  +
  xlab("Days") +
  scale_colour_manual(values = c("skyblue", "orange"))

# decline in effectiveness







