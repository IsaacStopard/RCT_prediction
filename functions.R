### calculating the top ups required to give the same mean
gen_cov <- function(retention, base, t = seq(0,24*30, 0.1)){
  data.frame("t" = t,
             "net_cov" = base*exp(-(1 - exp(-1/(retention*365)))*t))
}

# all top ups are given to people without nets

# calculate the mean coverage with declining bed net use
calc_mean_cov <- function(retention, base, min_t, max_t){
  rate <- 1 - exp(-1/(retention))
  int_u <- -1/rate * exp(-rate * max_t) # integral
  int_l <- -1/rate * exp(-rate * min_t) #
  mean <- (int_u - int_l) * base /(max_t - min_t)
  return(mean)
}

calc_base_mean <- function(retention, mean, min_t, max_t){
  rate <- 1 - exp(-1/(retention))
  int_u <- -1/rate * exp(-rate * max_t)
  int_l <- -1/rate * exp(-rate * min_t)
  base <- mean * (max_t - min_t) / (int_u - int_l) 
  return(base)
}

# decline in net efficacy - killing probability 
decline_d0 <- function(t, dn0, gamma_n){
  rate = -log(1/2)/gamma_n
  return(dn0 * exp(-t * rate))
}

# decline in net efficacy - repelling probability
decline_r0 <- function(t, rn0, rnm, gamma_n){
  rate = -log(1/2)/gamma_n
  return((rn0 - rnm) * exp(-t * rate) + rnm)
}

