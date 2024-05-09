rm(list = ls())

# find proportion of people with bednets
suppressPackageStartupMessages(library(ggplot2))
library(malariasimulation)
library(malariaEquilibrium)
library(reshape2)
library(tidyverse)

year <- 365
month <- 30
sim_length <- 3 * year
human_population <- 1000
starting_EIR <- 50

simparams <- get_parameters(
  list(
    human_population = human_population,
    model_seasonality = TRUE, # Let's try a bi-modal model
    g0 = 0.28605,
    g = c(0.20636, -0.0740318, -0.0009293),
    h = c(0.173743, -0.0730962, -0.116019),
    severe_incidence_rendering_min_ages = 2*year,
    severe_incidence_rendering_max_ages = 10*year
  )
)

simparams <- set_equilibrium(simparams, starting_EIR)

bednet_events = data.frame(
  timestep = c(1),
  name=c("Bednets 1")
)

retention <- 5 * year

# bednetparams <- set_bednets(
#   simparams,
#   timesteps = bednet_events$timestep,
#   coverages = c(1),
#   retention = c5 * year,
#   dn0 = matrix(c(.533), nrow=1, ncol=1),
#   rn = matrix(c(.56), nrow=1, ncol=1),
#   rnm = matrix(c(.24), nrow=1, ncol=1),
#   gamman = c(2.64 * 365))

bednetparams <- set_bednets(
  simparams,
  timesteps = bednet_events$timestep,
  coverages = c(1),
  retention = 5 * year,
  dn0 = matrix(c(.533), nrow=1, ncol=1),
  rn = matrix(c(.56), nrow=1, ncol=1),
  rnm = matrix(c(.24), nrow=1, ncol=1),
  gamman = rep(2.64 * 365, 1)
)

correlationsb1 <- get_correlation_parameters(bednetparams)
correlationsb1$inter_round_rho('bednets', 1) # if one nets are given to people who already have nets

output <- run_simulation(sim_length, bednetparams, correlations = correlationsb1)

plot(output$timestep, output$n_use_net/1000, type="l")

rate <- (1 - exp(-1/retention))

r_plot <- data.frame("t" = output$timestep,
           "n_use" = output$n_use_net/1000,
           "n_pred" = c(0,exp(-rate*(output$timestep-1))[-length(output$timestep)]))

ggplot(data = r_plot, aes(x = t)) +
  geom_line(aes(y = n_use), col = "blue", size = 1.5) +
  geom_line(aes(y = n_pred), size = 1.5) +
  ylab("% using nets") + xlab("Days")

r_fun <- function(retention){
  rate <- (1 - exp(-1/retention))
  return(data.frame("t" = seq(1, 365, 0.5)) %>% 
           mutate(r = exp(-rate*t),
                  retention = retention))
         
        }

r_df <- bind_rows(lapply(seq(1, 10, 0.5)*365, r_fun))

ggplot(data = r_df, aes(x = t, y = r, col = retention, group = factor(retention))) +
  geom_line(size = 1.5)

### adding new nets
bednet_events = data.frame(
  timestep = c(1, 2) * year,
  name=c("Bednets 1", "Bednets 2")
)

bednetparams <- set_bednets(
  bednetparams,
  timesteps = bednet_events$timestep,
  coverages = c(.8, .8),
  retention = c(5 * year, 1 * year),
  dn0 = matrix(c(.533, .45), nrow=2, ncol=1),
  rn = matrix(c(.56, .5), nrow=2, ncol=1),
  rnm = matrix(c(.24, .24), nrow=2, ncol=1),
  gamman = rep(2.64 * 365, 2)
)

correlationsb1 <- get_correlation_parameters(bednetparams)
correlationsb1$inter_round_rho('bednets', 0) # if one nets are given to people who already have nets

output <- run_simulation(sim_length, bednetparams, correlations = correlationsb1)


r_plot <- data.frame("t" = output$timestep,
                     "n_use" = output$n_use_net/1000)

ggplot(data = r_plot, aes(x = t, y = n_use)) +
  geom_line()

# setting multiple bed net events

# calculating the increase in nets with a random allocation (correlation parameters = 0)
r_plot$n_use[365*2]

0.642 + 0.8 - (0.8 * 0.642)

r_plot$n_use[365*2+1]

# trying different retention times
### adding new nets
bednet_events_1 = data.frame(
  timestep = c(1) * year,
  name=c("Bednets 1")
)

bednet_events_2 = data.frame(
  timestep = c(2) * year,
  name=c("Bednets 2")
)

bednetparams <- set_bednets(
  simparams,
  timesteps = bednet_events_1$timestep,
  coverages = c(.8),
  retention = 5 * year,
  dn0 = matrix(c(.533), nrow=1, ncol=1),
  rn = matrix(c(.56), nrow=1, ncol=1),
  rnm = matrix(c(.24), nrow=1, ncol=1),
  gamman = rep(2.64 * 365, 1)
)

bednetparams <- set_bednets(
  bednetparams,
  timesteps = bednet_events_2$timestep,
  coverages = c(.8),
  retention = 5 * year,
  dn0 = matrix(c(.533), nrow=1, ncol=1),
  rn = matrix(c(.56), nrow=1, ncol=1),
  rnm = matrix(c(.24), nrow=1, ncol=1),
  gamman = rep(2.64 * 365, 1)
)

correlationsb1 <- get_correlation_parameters(bednetparams)
correlationsb1$inter_round_rho('bednets', 0) # if one nets are given to people who already have nets

output <- run_simulation(sim_length, bednetparams, correlations = correlationsb1)

r_plot <- data.frame("t" = output$timestep,
                     "n_use" = output$n_use_net/1000)

ggplot(data = r_plot, aes(x = t, y = n_use)) +
  geom_line()





