library(invgamma)
library(MASS)
library(tidyverse)

generate_data <- function(n, p, m, signal = 10) {
  # n x p
  x <- mvrnorm(n, mu = rnorm(p),
               Sigma = diag(.1, nrow = p))
  
  x_tilde <-  mvrnorm(n, mu = rnorm(p),
                      Sigma = diag(.1, nrow = p))
  
  # n x 1
  
  eps <- rnorm(n, mean = 0, sd = rinvgamma(n, 100, 10))
  
  
  # p x 1
  beta <- c(rep(0, m), 
            rnorm((p - m) / 2, mean = signal, sd = .1),
            rnorm((p - m) / 2, mean = - signal, sd = .1)
  )
  y <- x %*% beta + eps
  
  return(list(
    x = x,
    x_tilde = x_tilde,
    y = y
  ))
}



# calculate the threshold in knock-off 
cal_t <- function(w, alpha) {
  for (t in sort(abs(w))){
    prop <- sum(w <= -t) / sum(w >= t)
    if (prop <= alpha) {
      break
    }
  }
  return(t)
}


cal_p <- function(p_vals, alpha, p) {
  # calculate the cutoff in BH
  bh_df <- data.frame(p_vals) %>%
    rownames_to_column(var = "original_index") %>%
    arrange(p_vals) %>%
    rownames_to_column(var = "index") %>%
    mutate(index = as.numeric(index),
           cutoff = p_vals - index / p * alpha)
  
  thres_bh <- bh_df %>%
    filter(cutoff <= 0) %>%
    dplyr::select(p_vals) %>%
    unlist %>% max
  return(thres_bh)
}




simulate <- function(n, p, m, K = 100, signal = 1) {
  fdps <- rep(NA, K)
  tdps <- rep(NA, K)
  fdps_bh <- rep(NA, K)
  tdps_bh <- rep(NA, K)
  for (k in 1:K) {
    set.seed(k)
    # generate data
    data_list <- generate_data(n, p, m, signal = signal)
    x <- data_list$x
    x_tilde <- data_list$x_tilde
    y <- data_list$y
    
    p_vals <- c()
    t_stats <- c()
    t_tilde_stats <- c()
    w_stats <- c()
    for (j in 1:p) {
      mod <- lm(y ~ x[, j] + x_tilde[, j])
      s = summary(mod)
      t_stats <- c(t_stats, s$coefficients[2, 3])
      p_vals <- c(p_vals, s$coefficients[2, 4])
      t_tilde_stats <- c(t_tilde_stats, s$coefficients[3, 3])
      
      w_stats <- c(w_stats, abs(t_stats[j]) - abs(t_tilde_stats[j] ))
    }
    
    thres_kf <- cal_t(w_stats, alpha)
    # false discovery proportion
    fdps[k] <- sum(w_stats[1:m] >= thres_kf) / p
    #print(fdp)
    
    # power
    tdps[k] <-sum(w_stats[(m + 1):p] >= thres_kf) / (p - m)
    #print(power)
    
    thres_bh <- cal_p(p_vals, alpha, p)
    fdps_bh[k] <- sum(p_vals[1:m] <= thres_bh) / p
    tdps_bh[k] <-sum(p_vals[(m + 1):p] <= thres_bh) / (p - m)
  }
  
  rates <- data.frame(
    fdp = c(fdps, fdps_bh),
    tdp = c(tdps, tdps_bh),
    method = rep(c("model-x knockoff", "BH"), each = K)
  )
  return(rates)
}

