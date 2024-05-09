rm(list = ls())

suppressPackageStartupMessages(library(ggplot2)); library(malariasimulation); library(malariaEquilibrium)
library(reshape2); library(tidyverse); library(readxl); library(rstan); library(pracma)
library(foresite); library(doParallel); library(foreach);
library(patchwork); library(pammtools); library(ggpattern);
library(cowplot)

source(file = "utils/functions.R"); source(file = "utils/retention_fit_top_up.R"); source(file = "utils/run_model_functions.R");
source(file = "read_data.R")

################################################
##### calculating the top up net coverages #####
################################################

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

saveRDS(top_up_all, file = "parameters/top_up_fits.rds")

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
