# rm(list = ls())
# 
# # find proportion of people with bednets
# suppressPackageStartupMessages(library(ggplot2)); library(malariasimulation); library(malariaEquilibrium)
# library(reshape2); library(tidyverse); library(readxl); library(rstan); library(pracma)
# 
# source(file = "functions.R")
# 
# ### RCT nets retention
# df <- read_excel("data/data_RCT_Benin.xlsx", sheet = "df") %>% 
#   mutate(Bed_net_use_both = Bed_net_use_both/100,
#          Bed_net_use_RCT = Bed_net_use_RCT/100)

binom_model <- "
data {
int<lower=0> N;
int n[N]; // number with net
vector[N] t; // time
int n_tested[N]; 

int<lower=0> pred_N;
vector[pred_N] pred_t;
}

parameters {
real<lower=0, upper = 40> retention;
real<lower = 0, upper=1> base; // initial number of people with RCT nets
}

transformed parameters {
real rate;
vector<lower=0, upper=1>[N] prob;

rate = (1 - exp(-1/(retention*365)));

prob = base * exp(-rate * t);
}

model {
// priors
retention ~ gamma(10, 10.0/2.5);
base ~ beta(1.5, 1.5);

// likelihood
n ~ binomial(n_tested, prob);   
}

generated quantities {
vector[pred_N] pred_p;

pred_p = base * exp(-rate * pred_t);
}
"

binom_model_base <- "
data {
int<lower=0> N;
int n[N]; // number with net
vector[N] t; // time
int n_tested[N];

int<lower=0> pred_N;
vector[pred_N] pred_t;
real<lower = 0, upper=1> base; // initial number of people with RCT nets
}

parameters {
real<lower=0, upper = 40> retention;
}

transformed parameters {
real rate;
vector<lower=0, upper=1>[N] prob;

rate = (1 - exp(-1/(retention*365)));

prob = base * exp(-rate * t);
}

model {
// priors
retention ~ gamma(10, 10.0/2.5);

// likelihood
n ~ binomial(n_tested, prob);   
}

generated quantities {
vector[pred_N] pred_p;

pred_p = base * exp(-rate * pred_t);
}
"

# model of the decline in all net types


multinom_model <- "
data {
  int K; // number of outcomes (3)
  int N; // number of times
  vector[N] t; // time
  int<lower=0> n[N, K]; // number with (1) rct net, (2) standard net & (3) no net
  
  int<lower=0> pred_N; // for generated quantities
  vector[pred_N] pred_t;
  
  //real<lower = 0, upper=1> base; // initial number of people with RCT nets
}

parameters {
  real<lower=0, upper = 40> retention;
  simplex[K] base; // initial proportion with (1) rct net, (2) standard net & (3) no net
}

transformed parameters {
  simplex[K] prob[N];
  
  real rate;
  rate = (1 - exp(-1/(retention*365)));
  
  for(i in 1:N){
    prob[i] = [base[1] * exp(-rate * t[i]), 
               base[2] * exp(-rate * t[i]),
               1.0 - base[1] * exp(-rate * t[i]) - base[2] * exp(-rate * t[i])]';
    }
}

model {
  // priors
  retention ~ gamma(10, 10.0/2.5);
  base ~ dirichlet(rep_vector(1.0, K));
  
  for(i in 1:N){
    target += multinomial_lpmf(n[i,] | prob[i]);
  }
 }

generated quantities {
  vector[pred_N] pred_p_rct;
  vector[pred_N] pred_p_p_only;
  vector[pred_N] pred_p_both;
  
  pred_p_rct = base[1] * exp(-rate * pred_t);
  pred_p_p_only = base[2] * exp(-rate * pred_t);
  pred_p_both = pred_p_rct + pred_p_p_only;

}
"

# compiling the models
stanDso_in <- rstan::stan_model(model_code = binom_model) 

stanDso_base_in <- rstan::stan_model(model_code = binom_model_base) 

stanDso_base_both_in <- rstan::stan_model(model_code = multinom_model) 

# data

get_top_up_multi <- function(net,
                             country,
                             df,
                             stanDso_base_both = stanDso_base_both_in){
  
  # subsetting the data so it's only for a specific arm
  # putting the data into long format depending on the net type
  subset_data <- subset(df, Net == net & Country == country)
  
  c_data <-  subset_data %>% 
    gather(key = "type", value = "Bed_net_use", Bed_net_use_both, Bed_net_use_RCT) %>% 
    mutate(net = ifelse(type == "Bed_net_use_both",1,2))
  
  # getting the rct data only
  rct_data <- subset(c_data, type == "Bed_net_use_RCT" & Time_months != 0)
  
  # calculating the mean bed net use of both nets
  mean_bed_net_use_both <-  weighted.mean(x = subset_data$Bed_net_use_both, w = subset_data$Bed_net_n_tested, na.rm = TRUE)
  
  pred_t <- seq(0, max(c_data$Time_months*30), 0.5)
  
  data_in <- list(K = 3,
                  N = nrow(rct_data),
                  n = as.matrix(data.frame(rct_data$Bed_net_n_use_RCT,
                                           rct_data$Bed_net_n_use_both - rct_data$Bed_net_n_use_RCT,
                                           rct_data$Bed_net_n_tested - rct_data$Bed_net_n_use_both)),
                  t = rct_data$Time_months*30, # assuming 30 days in a month
                  pred_N = length(pred_t),
                  pred_t = pred_t
  )
  
  finit <- function(){
    list(retention = 5,
         base = c(0.5, 0.3, 0.2)
    )
  }
  
  fit <- rstan::sampling(stanDso_base_both, 
                         data = data_in, 
                         iter = 2500, warmup=1250,
                         init= finit, cores = 4) 
  
  fit_e <- rstan::extract(fit)
  plot_df <- data.frame("t" = pred_t,
                        "pred_p" = apply(fit_e$pred_p_rct, 2, quantile, prob = 0.5),
                        "pred_p_l" = apply(fit_e$pred_p_rct, 2, quantile, prob = 0.025),
                        "pred_p_u" = apply(fit_e$pred_p_rct, 2, quantile, prob = 0.975),
                        "pred_p_b" = apply(fit_e$pred_p_both, 2, quantile, prob = 0.5),
                        "pred_p_b_l" = apply(fit_e$pred_p_both, 2, quantile, prob = 0.025),
                        "pred_p_b_u" = apply(fit_e$pred_p_both, 2, quantile, prob = 0.975))
  
  plot <- ggplot(data = plot_df) +
    geom_ribbon(aes(x = t, ymin = pred_p_l, ymax = pred_p_l), fill = "skyblue", alpha = 0.75) +
    geom_ribbon(aes(x = t, ymin = pred_p_b_l, ymax = pred_p_b_l), alpha = 0.75) +
    geom_line(aes(x = t, y = pred_p), col = "skyblue", linewidth = 0.5) +
    geom_line(aes(x = t, y = pred_p_b), linewidth = 0.5) +
    geom_point(data = rct_data, aes(x = Time_months * 30, y = Bed_net_n_use_RCT/Bed_net_n_tested), col = "skyblue", size = 2.5) +
    geom_point(data = rct_data, aes(x = Time_months * 30, y = Bed_net_n_use_both/Bed_net_n_tested), size = 2.5) +
    theme_bw() +
    ylim(c(0, 1))
  
  return(list("rct_data" = rct_data, "mean_bed_net_use_both" = mean_bed_net_use_both, "fit" = fit,
              "plot" = plot))
  
}


get_top_up <- function(net,
                       country,
                       df,
                       tu_diff_time, # top ups every 6 months
                       sim_length,
                       int_bed_net_time,
                       int_date,
                       stanDso = stanDso_in,
                       stanDso_base = stanDso_base_in){
  
  set.seed(12345)
  
  if(tu_diff_time < 2){stop("tu_diff_time must be greater than 1")}
  # subsetting the data so it's only for a specific arm
  # putting the data into long format depending on the net type
  subset_data <- subset(df, Net == net & Country == country)
  subset_data$t <- as.double(difftime(subset_data$date, int_date))
  
  
  c_data <-  subset_data %>% 
  gather(key = "type", value = "Bed_net_use", Bed_net_use_both, Bed_net_use_RCT) %>% 
    mutate(net = ifelse(type == "Bed_net_use_both",1,2))
  
  # getting the rct data only
  rct_data <- subset(c_data, type == "Bed_net_use_RCT" & is.na(Time_months) != 1 & Time_months != 30) # excluding the 30 months
  
  
  # calculating the mean bed net use of both nets
  mean_bed_net_use_both <-  weighted.mean(x = subset(subset_data, is.na(Time_months) != 1)$Bed_net_use_both, 
                                          w = subset(subset_data, is.na(Time_months) != 1)$Bed_net_n_tested, 
                                          na.rm = TRUE)
  
  pred_t <- seq(0, max(c_data$t), 0.5)
  
  data_in <- list(N = nrow(rct_data),
                n = rct_data$Bed_net_n_use_RCT,
                t = rct_data$t, # assuming 30 days in a month
                pred_N = length(pred_t),
                pred_t = pred_t,
                n_tested = rct_data$Bed_net_n_tested
                )
  
  finit <- function(){
  list(retention = 5,
       base = 0.5
  )
    }

  # fit <- rstan::sampling(stanDso, 
  #                        data = data_in, 
  #                        iter = 5000, warmup=2500,
  #                        init= finit, cores = 4) 
  
  # this model assumes the initial coverage is the same as the mean net use of both net types
  fit_base <- rstan::sampling(stanDso_base, 
                              data = append(data_in, list("base" = mean_bed_net_use_both)), 
                              iter = 5000, warmup=2500,
                              init= finit, cores = 4)

  #fit_e <- rstan::extract(fit)
  
  fit_e_base <- rstan::extract(fit_base)
  
  # extracting the predicted values
  
  plot_df_base <- data.frame("t" = pred_t,
                      "pred_p" = apply(fit_e_base$pred_p, 2, quantile, prob = 0.5),
                      "pred_p_l" = apply(fit_e_base$pred_p, 2, quantile, prob = 0.025),
                      "pred_p_u" = apply(fit_e_base$pred_p, 2, quantile, prob = 0.975))
  
  c_data$lower <- c_data$Bed_net_use - 1.96 * sqrt((c_data$Bed_net_use * (1 - c_data$Bed_net_use))/c_data$Bed_net_n_tested)
  c_data$upper <- c_data$Bed_net_use + 1.96 * sqrt((c_data$Bed_net_use * (1 - c_data$Bed_net_use))/c_data$Bed_net_n_tested)
  
  # give everyone pyrethroid nets
  # the day after give the trial nets replacing the pyrethroid net
  # calculate the initial pyrethroid net coverage as being the value that gives the same mean
  
  ###
  # initially coverage is 100%
  # finding the time until the first top up to give the same mean
  
  # time until initial top up
  
  # so the total nets don't decline but stays fixed at the mean value
  # calculates the proportion of nets that need to be added to observe a decline in the trial nets
  # assuming that all nets are given to people who already have nets
  covs_df <- data.frame("t" = seq(0, sim_length - int_bed_net_time + tu_diff_time, tu_diff_time)
                        ) %>% 
    mutate(int_bed_net_cov_m = NA,
           int_bed_net_cov_l = NA,
           int_bed_net_cov_u = NA,
           mean_cov_m = NA,
           mean_cov_l = NA,
           mean_cov_u = NA)
  
  # replacing all the nets with intervention nets
  covs_df[1,] <- c(0, rep(mean_bed_net_use_both, 3), rep(NA, 3))
  
  # covs_df <- covs_df %>% mutate(diff_int_mean_cov_m = c(0, rep(NA, nrow(covs_df)-1)),
  #                               diff_int_mean_cov_l = c(0, rep(NA, nrow(covs_df)-1)),
  #                               diff_int_mean_cov_u = c(0, rep(NA, nrow(covs_df)-1)))
  
  # calculating the differences between the nets lost at different times
  for(i in 2:nrow(covs_df)){
    covs_df[i, "int_bed_net_cov_m"] <- gen_cov(retention = median(fit_e_base$retention)*365,
                                             base = covs_df[1, "int_bed_net_cov_m"],
                                             t = covs_df[i, "t"])$net_cov
    
    # if lower then the longer retention is used
    covs_df[i, "int_bed_net_cov_l"] <- gen_cov(retention = quantile(fit_e_base$retention, probs = c(0.975))[[1]],
                                               base = covs_df[1, "int_bed_net_cov_l"],
                                               t = covs_df[i, "t"])$net_cov
    
    covs_df[i, "int_bed_net_cov_u"] <- gen_cov(retention = quantile(fit_e_base$retention, probs = c(0.025))[[1]],
                                               base = covs_df[1, "int_bed_net_cov_u"],
                                               t = covs_df[i, "t"])$net_cov
  }
  
  for(i in 1:(nrow(covs_df)-1)){
    covs_df[i, "mean_cov_m"] <- calc_mean_cov(retention = median(fit_e_base$retention),
                                              base = covs_df[1, "int_bed_net_cov_m"],
                                              min_t = covs_df[i, "t"],
                                              max_t = covs_df[(i + 1), "t"])
    
    covs_df[i, "mean_cov_l"] <- calc_mean_cov(retention = quantile(fit_e_base$retention, probs = c(0.025))[[1]],
                                              base = covs_df[1, "int_bed_net_cov_m"],
                                              min_t = covs_df[i, "t"],
                                              max_t = covs_df[(i + 1), "t"])
    
    covs_df[i, "mean_cov_u"] <- calc_mean_cov(retention = quantile(fit_e_base$retention, probs = c(0.975))[[1]],
                                              base = covs_df[1, "int_bed_net_cov_m"],
                                              min_t = covs_df[i, "t"],
                                              max_t = covs_df[(i + 1), "t"])
  }
  
  prop_t <- (tu_diff_time - 1)/tu_diff_time
  covs_df <- covs_df %>% mutate(top_up = mean_cov_m/prop_t,
                                top_up_l = mean_cov_l/prop_t,
                                top_up_u = mean_cov_u/prop_t)
  
  # covs_df <- covs_df[-nrow(covs_df),]
  # 
  # for(i in 2:nrow(covs_df)){
  #   covs_df[i, "diff_int_mean_cov_m"] = covs_df[(i-1), "mean_cov_m"] - covs_df[i, "mean_cov_m"]
  #   
  #   covs_df[i, "diff_int_mean_cov_l"] = covs_df[(i-1), "mean_cov_l"] - covs_df[i, "mean_cov_l"]
  #   
  #   covs_df[i, "diff_int_mean_cov_u"] = covs_df[(i-1), "mean_cov_u"] - covs_df[i, "mean_cov_u"]
  # }
  # 
  # covs_df <- covs_df %>% mutate(prop_int_m = mean_cov_m/covs_df[1, "mean_cov_m"],
  #                               prop_int_l = mean_cov_l/covs_df[1, "mean_cov_l"],
  #                               prop_int_u = mean_cov_u/covs_df[1, "mean_cov_u"],
  #                               top_up = c(0, rep(NA, nrow(covs_df)-1)),
  #                               top_up_l = c(0, rep(NA, nrow(covs_df)-1)),
  #                               top_up_u = c(0, rep(NA, nrow(covs_df)-1)))
  # 
  # for(i in 2:nrow(covs_df)){
  #   covs_df[i, "top_up"] <- covs_df[i,"diff_int_mean_cov_m"] / covs_df[(i-1), "prop_int_m"]
  #   covs_df[i, "top_up_l"] <- covs_df[i,"diff_int_mean_cov_l"] / covs_df[(i-1), "prop_int_l"]
  #   covs_df[i, "top_up_u"] <- covs_df[i,"diff_int_mean_cov_u"] / covs_df[(i-1), "prop_int_u"]
  #   
  # }
  
  plot <- ggplot(data = c_data[-which(is.na(c_data$Bed_net_use)),]) +
    #geom_ribbon(data = plot_df, aes(x = t, ymin = pred_p_l, ymax = pred_p_u), alpha = 0.25) +
    geom_ribbon(data = plot_df_base, aes(x = t, ymin = pred_p_l, ymax = pred_p_u), alpha = 0.25, fill = "skyblue") +
    geom_pointrange(aes(x = t,
                        y = Bed_net_use, ymin = lower, ymax = upper, fill = type), size = 0.8,
                    shape = 21) +
    scale_fill_manual(values = c("black", "skyblue"), labels = c("Any type of bed net", "Trial bed net"), name = "") +
    #geom_line(data = plot_df, aes(x = t, y = pred_p), linewidth = 1.5) +
    pammtools::geom_stepribbon(data = na.omit(covs_df), aes(x = t, ymin = mean_cov_l, ymax = mean_cov_u), alpha = 0.1) +
    geom_line(data = plot_df_base, aes(x = t, y = pred_p), linewidth = 1.1, col = "skyblue") +
    geom_step(data = na.omit(covs_df), aes(x = t, y = mean_cov_m), alpha = 0.5, linewidth = 0.5) +
    theme_bw() + theme(text = element_text(size = 18)) +
    xlab("Days since start of trial") +
    ylab("Proportion using bed net") +
    geom_hline(yintercept = mean_bed_net_use_both, linetype = 2, linewidth = 1) +
    scale_x_continuous(breaks = seq(-200, 1400, 200)) +
    coord_cartesian(xlim = c(-200, 1500)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1))
  
  return(list("rct_data" = rct_data, "mean_bed_net_use_both" = mean_bed_net_use_both, 
               "fit_base" = fit_base, #"fit" = fit,
                          "plot" = plot, "covs_df" = covs_df))
}
