rm(list = ls())

# find proportion of people with bednets
suppressPackageStartupMessages(library(ggplot2)); library(malariasimulation); library(malariaEquilibrium)
library(reshape2); library(tidyverse); library(readxl); library(rstan)

df <- read_excel("data/data_RCT_Benin.xlsx", sheet = "df") %>% 
  mutate(Bed_net_use_both = Bed_net_use_both/100,
         Bed_net_use_RCT = Bed_net_use_RCT/100)

binom_model <- "
data {
int<lower=0> N;
int n[N]; // number with net
vector[N] t; // time
int net[N]; // net type

int<lower=0> pred_N;
vector[pred_N] pred_t;
}

parameters {
real<lower=0, upper = 40> retention;
real<lower = 0, upper=1> base_c;
real<lower = 0, upper=1> base_a;
}

transformed parameters {
real rate;
vector<lower=0, upper=1>[N] prob;

rate = (1 - exp(-1/(retention*365)));

for(i in 1:N){
if(net[i] == 1){ // 1 is for both nets
  prob[i] = base_c * exp(-rate * t[i]) + base_a * exp(-rate * t[i]);
  
} else{
  prob[i] = base_a * exp(-rate * t[i]);
}


}
}

model {
// priors
retention ~ gamma(10, 10.0/2.5);
base_c ~ beta(1.5, 1.5);
base_a ~ beta(1.5, 1.5);

// likelihood
n ~ binomial(1372, prob);   
}

generated quantities {
vector[pred_N] pred_pc;
vector[pred_N] pred_pa;

pred_pc = base_c * exp(-rate * pred_t) + base_a * exp(-rate * pred_t);
pred_pa = base_a * exp(-rate * pred_t);
}
"

stanDso <- stan_model(model_code = binom_model) 
c_data <- subset(df, Location == "Cove") %>% 
  gather(key = "type", value = "Bed_net_use", Bed_net_use_both, Bed_net_use_RCT) %>% 
  mutate(net = ifelse(type == "Bed_net_use_both",1,2))

pred_t <- seq(0, max(c_data$Time_months*30), 0.5)

data_in <- list(N = nrow(c_data),
  n = round(c_data$Bed_net_use*1372, digits = 0),
  t = c_data$Time_months*30,
  net = as.integer(c_data$net),
  pred_N = length(pred_t),
  pred_t = pred_t
)

finit <- function(){
  list(retention = 5,
       c = 0.2,
       base = 0.5
  )
}

fit <- sampling(stanDso, data = data_in, iter = 5000, warmup=2500,
                init= finit) 

fit_e <- rstan::extract(fit)
mean(fit_e$retention)

plot_df <- data.frame("t" = pred_t,
           "pred_c" = apply(fit_e$pred_pc, 2, quantile, prob = 0.5),
           "pred_c_l" = apply(fit_e$pred_pc, 2, quantile, prob = 0.025),
           "pred_c_u" = apply(fit_e$pred_pc, 2, quantile, prob = 0.975),
           "pred_a" = apply(fit_e$pred_pa, 2, quantile, prob = 0.5),
           "pred_a_l" = apply(fit_e$pred_pa, 2, quantile, prob = 0.025),
           "pred_a_u" = apply(fit_e$pred_pa, 2, quantile, prob = 0.975))

gen_cov <- function(retention, base_a, base_c, t = seq(0,24*30, 0.1)){
  data.frame("t" = t,
             "tot_net_cov" = base_c*exp(-(1 - exp(-1/(retention*365)))*t),
             "rct_net_cov" = base_a*exp(-(1 - exp(-1/(retention*365)))*t))
}



c_df <- gen_cov(retention = mean(fit_e$retention),
        base_c = mean(fit_e$base_c),
        base_a = mean(fit_e$base_a))



ggplot(data = c_data %>% mutate(prob = apply(fit_e$prob, 2, mean))) +
  geom_ribbon(data = plot_df, aes(x = t, ymin = pred_c_l, ymax = pred_c_u), alpha = 0.25, fill = "blue") +
  geom_ribbon(data = plot_df, aes(x = t, ymin = pred_a_l, ymax = pred_a_u), alpha = 0.25, fill = "grey50") +
  
  geom_point(aes(x = Time_months*30,
                 y = Bed_net_use, 
                 col = factor(net)), size = 3) +
  
  scale_colour_manual(name = "Nets", values = c("blue", "grey50"), labels = c("all", "RCT nets")) +
  
  geom_line(data = plot_df,
            aes(x = t, y = pred_c), col = "blue", size = 1.5) +
  geom_line(data = plot_df,
            aes(x = t, y = pred_a), col = "grey50", size = 1.5) +
  theme_bw() +
  xlab("Days since start of trial") +
  ylab("Proportion using bed net") +
  ylim(0.4, 1)



  ggplot(data = df %>% gather(key = "Pop", value = "Bed_net_use", Bed_net_use_both, Bed_net_use_RCT), 
       aes(x = Time_months, y = Bed_net_use, col = Pop)) +
  geom_point(size = 3) +
  facet_wrap(~Location)
  