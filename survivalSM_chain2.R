#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(invgamma)
library(dplyr)
library(mvtnorm)
library(stats4)
set.seed(12345 + as.numeric(args[1]))


df1 <- read.csv(paste0("simudata", as.numeric(args[1]), ".csv"))

niter <- 50000

lag = 1
N = nrow(df1)
Ni = length(unique(df1$pid))

y = df1$FEV1

ni <- aggregate(df1$lamdxage, by = list(df1$pid), length)$x

df <- df1[c(2, 3, 8, 5, 9, 10, 12, 13)]

ni_l = c(1, (cumsum(ni)+1)[-length(ni)])
ni_u = cumsum(ni)
s_df <- df %>% group_by(pid) %>% slice(which.max(timeSince_t0))
pt = unique(df$pid)

df_t_Ti <- df[1:2]
colnames(df_t_Ti) <- c("pid", "t_Ti")

gamma_h01 = c()
gamma_h02 = c()

beta1 = c()

gam1 = c()


main_para <- readRDS(paste0("longitudinalSM_mainpara_estimates_simulation", as.numeric(args[1]),".rds"))

alpha0 =  main_para["alpha0", "Mean"]
alpha1 = main_para["alpha1", "Mean"]

sigmasq = main_para["sigmasq", "Mean"]

tausq = main_para["tausq", "Mean"]

re_stoc <- readRDS(paste0("longitudinalSM_random_terms_estimates_simulation",as.numeric(args[1]),".rds"))

Wi <- re_stoc[c(1:N), "Mean"]
ui <- re_stoc[-c(1:N), "Mean"]

H <- main_para["H", "Mean"]
kappa <- main_para["kappa", "Mean"]

gamma_h01[1] = as.numeric(args[8])
gamma_h02[1] = as.numeric(args[9])

beta1[1] = as.numeric(args[10])

gam1[1] = as.numeric(args[11])


#MH for gamma_h01
f_gamma_h01 <- function(ui, alpha0, alpha1, gamma_h01, gamma_h02, beta1, gam1, Wi, df_t_Ti, y, ni, N, X, W, deltai, sd_gamma_h01_star)
{
  p_gamma_h01 <- function(val)
  {
    df_t_Ti$diff_t <- df_t_Ti %>% group_by(pid) %>% mutate(diff_t_Ti = exp(beta1*alpha1*t_Ti) - exp(beta1*alpha1*lag(t_Ti))) %>% pull(diff_t_Ti)
    
    df_t_Ti$Wi <- Wi
    
    df_t_Ti <- df_t_Ti[complete.cases(df_t_Ti),]
    
    l <- df_t_Ti %>% group_by(pid) %>% mutate(ll =  t_Ti^(val - 1)*exp(beta1*Wi)*diff_t/(beta1*alpha1)) %>% summarise(l = sum(ll), .groups = 'drop') %>% pull(l)
    
    return(sum(deltai*log(val) - deltai*val*log(gamma_h02) + deltai*(val - 1)*log(W[,5]) - val*gamma_h02^(-val)*exp(gam1*W[,1] + beta1*(ui + alpha0))*l) + dgamma(val, shape = 0.01, rate = 0.01, log = T))
  }
  
  p_gamma_h01_mle <- function(val)
  {
    df_t_Ti$diff_t <- df_t_Ti %>% group_by(pid) %>% mutate(diff_t_Ti = exp(beta1*alpha1*t_Ti) - exp(beta1*alpha1*lag(t_Ti))) %>% pull(diff_t_Ti)
    
    df_t_Ti$Wi <- Wi
    
    df_t_Ti <- df_t_Ti[complete.cases(df_t_Ti),]
    
    l <- df_t_Ti %>% group_by(pid) %>% mutate(ll =  t_Ti^(val - 1)*exp(beta1*Wi)*diff_t/(beta1*alpha1)) %>% summarise(l = sum(ll), .groups = 'drop') %>% pull(l)
    
    return(-(sum(deltai*log(val) - deltai*val*log(gamma_h02) + deltai*(val - 1)*log(W[,5]) - val*gamma_h02^(-val)*exp(gam1*W[,1] + beta1*(ui + alpha0))*l) + dgamma(val, shape = 0.01, rate = 0.01, log = T)))
  }
  
  #gamma_h01_post <- stats4::mle(p_gamma_h01_mle, start = list(val = gamma_h01), method = "Brent", lower = 1e-11, upper = 10)@coef
  
  #gamma_h01_star <- rlnorm(1, meanlog = log(gamma_h01_post), sdlog = sd_gamma_h01_star)
  
  gamma_h01_star <- rlnorm(1, meanlog = log(gamma_h01), sdlog = sd_gamma_h01_star)
  
  u_gamma_h01 <- runif(1)
  
  #log_r_gamma_h01 <- p_gamma_h01(gamma_h01_star) + dlnorm(gamma_h01, meanlog = log(gamma_h01_post), sdlog = sd_gamma_h01_star, log = T) - p_gamma_h01(gamma_h01) - dlnorm(gamma_h01_star, meanlog = log(gamma_h01_post), sdlog = sd_gamma_h01_star, log = T)
  
  log_r_gamma_h01 <- p_gamma_h01(gamma_h01_star) + dlnorm(gamma_h01, meanlog = log(gamma_h01_star), sdlog = sd_gamma_h01_star, log = T) - p_gamma_h01(gamma_h01) - dlnorm(gamma_h01_star, meanlog = log(gamma_h01), sdlog = sd_gamma_h01_star, log = T)
  
  if(min(exp(log_r_gamma_h01), 1) >= u_gamma_h01){
    update_gamma_h01 <- 1
    return(list(gamma_h01_star, update_gamma_h01))
  } else 
  {
    update_gamma_h01 <- 0
    return(list(gamma_h01, update_gamma_h01))
  }
}

update_gamma_h01 <- c()

#for(i in 2:1000)
#{
#  res_gamma_h01 <- f_gamma_h01(ui = ui, alpha0, alpha1, gamma_h01 = gamma_h01[i-1], gamma_h02 = gamma_h02[1], beta1, gam1, Wi = Wi, df_t_Ti, y = df$FEV1, ni, N, X = df[1:6], W = s_df[c(3:6, 2)], deltai = s_df$event, sd_gamma_h01_star = 0.05)

#  gamma_h01[i] <- res_gamma_h01[[1]]

#  update_gamma_h01[i] <- res_gamma_h01[[2]]
#}

#MH for gamma_h02
f_gamma_h02 <- function(ui, alpha0, alpha1, gamma_h01, gamma_h02, beta1, gam1, Wi, df_t_Ti, y, ni, N, X, W, deltai, sd_gamma_h02_star)
{
  p_gamma_h02 <- function(val)
  {
    df_t_Ti$diff_t <- df_t_Ti %>% group_by(pid) %>% mutate(diff_t_Ti = exp(beta1*alpha1*t_Ti) - exp(beta1*alpha1*lag(t_Ti))) %>% pull(diff_t_Ti)
    
    df_t_Ti$Wi <- Wi
    
    df_t_Ti <- df_t_Ti[complete.cases(df_t_Ti),]
    
    l <- df_t_Ti %>% group_by(pid) %>% mutate(ll =  t_Ti^(gamma_h01 - 1)*exp(beta1*Wi)*diff_t/(beta1*alpha1)) %>% summarise(l = sum(ll), .groups = 'drop') %>% pull(l)
    
    return(sum(-gamma_h01*deltai*log(val) - gamma_h01*val^(-gamma_h01)*exp(gam1*W[,1] + beta1*(ui + alpha0))*l) + dgamma(val, shape = 0.1, rate = 0.1, log = T))
  }
  
  p_gamma_h02_mle <- function(val)
  {
    df_t_Ti$diff_t <- df_t_Ti %>% group_by(pid) %>% mutate(diff_t_Ti = exp(beta1*alpha1*t_Ti) - exp(beta1*alpha1*lag(t_Ti))) %>% pull(diff_t_Ti)
    
    df_t_Ti$Wi <- Wi
    
    df_t_Ti <- df_t_Ti[complete.cases(df_t_Ti),]
    
    l <- df_t_Ti %>% group_by(pid) %>% mutate(ll =  t_Ti^(gamma_h01 - 1)*exp(beta1*Wi)*diff_t/(beta1*alpha1)) %>% summarise(l = sum(ll), .groups = 'drop') %>% pull(l)
    
    return(-(sum(-gamma_h01*deltai*log(val) - gamma_h01*val^(-gamma_h01)*exp(gam1*W[,1] + beta1*(ui + alpha0))*l) + dgamma(val, shape = 0.01, rate = 0.01, log = T)))
  }
  #gamma_h02_post <- stats4::mle(p_gamma_h02_mle, start = list(val = gamma_h02), method = "Brent", lower = 1e-11, upper = 10)@coef
  
  #gamma_h02_star <- rlnorm(1, meanlog = log(gamma_h02_post), sdlog = sd_gamma_h02_star)
  gamma_h02_star <- rlnorm(1, meanlog = log(gamma_h02), sdlog = sd_gamma_h02_star)
  u_gamma_h02 <- runif(1)
  
  #log_r_gamma_h02 <- p_gamma_h02(gamma_h02_star) + dlnorm(gamma_h02, meanlog = log(gamma_h02_post), sdlog = sd_gamma_h02_star, log = T) - p_gamma_h02(gamma_h02) - dlnorm(gamma_h02_star, meanlog = log(gamma_h02_post), sdlog = sd_gamma_h02_star, log = T)
  
  log_r_gamma_h02 <- p_gamma_h02(gamma_h02_star) + dlnorm(gamma_h02, meanlog = log(gamma_h02_star), sdlog = sd_gamma_h02_star, log = T) - p_gamma_h02(gamma_h02) - dlnorm(gamma_h02_star, meanlog = log(gamma_h02), sdlog = sd_gamma_h02_star, log = T)
  
  r_gamma_h02 = exp(log_r_gamma_h02)
  if(min(r_gamma_h02, 1) >= u_gamma_h02){
    update_gamma_h02 <- 1
    return(list(gamma_h02_star, update_gamma_h02))
  } else 
  {
    update_gamma_h02 <- 0
    return(list(gamma_h02, update_gamma_h02))
  }
}

update_gamma_h02 <- c()

#for(i in 2:1000)
#{
#  res_gamma_h02 <- f_gamma_h02(ui = ui, alpha0, alpha1, gamma_h01 = gamma_h01[1], gamma_h02_current = gamma_h02[i-1], beta1, gam1, Wi = Wi, df_t_Ti, y = df$FEV1, ni, N, X = df[1:6], W = s_df[c(3:6, 2)], deltai = s_df$event, sd_gamma_h02_star = 1)

# gamma_h02[i] <- res_gamma_h02[[1]]

# update_gamma_h02[i] <- res_gamma_h02[[2]]
#}

#MH for beta1
f_beta1 <- function(ui, alpha0, alpha1, gamma_h01, gamma_h02, beta1, gam1, Wi, df_t_Ti, y, ni, N, X, W, deltai, sd_beta1_star)
{
  p_beta1 <- function(val)
  {
    df_t_Ti$diff_t <- df_t_Ti %>% group_by(pid) %>% mutate(diff_t_Ti = exp(val*alpha1*t_Ti) - exp(val*alpha1*lag(t_Ti))) %>% pull(diff_t_Ti)
    
    df_t_Ti$Wi <- Wi
    df_t_Ti <- df_t_Ti[complete.cases(df_t_Ti),]
    
    l <- df_t_Ti %>% group_by(pid) %>% mutate(ll =  t_Ti^(gamma_h01-1)*exp(val*Wi)*diff_t/(val*alpha1)) %>% summarise(l = sum(ll), .groups = 'drop') %>% pull(l)
    
    return(sum(deltai*val*(ui + alpha0 + alpha1*W[,5] + Wi[ni_u]) - gamma_h01*gamma_h02^(-gamma_h01)*exp(gam1*W[,1] + val*(ui + alpha0))*l) + dnorm(val, mean = 0, sd = 100, log = T))
  }
  
  p_beta1_mle <- function(val)
  {
    df_t_Ti$diff_t <- df_t_Ti %>% group_by(pid) %>% mutate(diff_t_Ti = exp(val*alpha1*t_Ti) - exp(val*alpha1*lag(t_Ti))) %>% pull(diff_t_Ti)
    
    df_t_Ti$Wi <- Wi
    df_t_Ti <- df_t_Ti[complete.cases(df_t_Ti),]
    
    l <- df_t_Ti %>% group_by(pid) %>% mutate(ll =  t_Ti^(gamma_h01-1)*exp(val*Wi)*diff_t/(val*alpha1)) %>% summarise(l = sum(ll), .groups = 'drop') %>% pull(l)
    
    return(-(sum(deltai*val*(ui + alpha0 + alpha1*W[,5] + Wi[ni_u]) - gamma_h01*gamma_h02^(-gamma_h01)*exp(gam1*W[,1] + val*(ui + alpha0))*l) + dnorm(val, mean = 0, sd = 100, log = T)))
  }
  
  #beta1_post <- stats4::mle(p_beta1_mle, start = list(val = beta1), method = "Brent", lower = -4, upper = 4)@coef
  
  #beta1_star <- rnorm(1, mean = beta1_post, sd = sd_beta1_star)
  
  beta1_star <- rnorm(1, mean = beta1, sd = sd_beta1_star)
  u_beta1 <- runif(1)
  
  log_r_beta1 <- p_beta1(beta1_star) + dnorm(beta1, mean = beta1_star, sd = sd_beta1_star, log = T) - p_beta1(beta1) - dnorm(beta1_star, mean = beta1, sd = sd_beta1_star, log = T)
  
  if(min(exp(log_r_beta1), 1) >= u_beta1)
  {
    update_beta1 <- 1
    return(list(beta1_star, update_beta1))
  } else {
    update_beta1 <- 0
    return(list(beta1, update_beta1))
  }
}
update_beta1 <- c()

#for(i in 2:1000)
#{
#  res_beta1 <- f_beta1(ui = ui[,1], alpha0, alpha1, gamma_h01, gamma_h02, beta1 = beta1[i - 1], beta2 = beta2, gam1, Wi = Wi[,1], Bi = Bi[,1], df_t_Ti, y = df$FEV1, ni, N, X = df[1:6], W = s_df[c(3:6, 2)], deltai = s_df$event, sd_beta1_star = 0.1)

#  beta1[i] <- res_beta1[[1]]
#  update_beta1[i] <- res_beta1[[2]]
#}

#MH for gamma1
f_gam1 <- function(ui, alpha0, alpha1, gamma_h01, gamma_h02, beta1, gam1, Wi, df_t_Ti, y, ni, N, X, W, deltai, sd_gam1_star)
{
  
  p_gam1 <- function(val)
  {
    df_t_Ti$diff_t <- df_t_Ti %>% group_by(pid) %>% mutate(diff_t_Ti = exp(beta1*alpha1*t_Ti) - exp(beta1*alpha1*lag(t_Ti))) %>% pull(diff_t_Ti)
    
    df_t_Ti$Wi <- Wi
    
    df_t_Ti <- df_t_Ti[complete.cases(df_t_Ti),]
    
    l <- df_t_Ti %>% group_by(pid) %>% mutate(ll =  t_Ti^(gamma_h01-1)*exp(beta1*Wi)*diff_t/(beta1*alpha1)) %>% summarise(l = sum(ll), .groups = 'drop') %>% pull(l)
    
    return(sum(deltai*val*W[,1] - gamma_h01*gamma_h02^(-gamma_h01)*exp(val*W[,1] + beta1*(ui + alpha0))*l) + dnorm(val, mean = 0, sd = 100, log = T))
  }
  
  p_gam1_mle <- function(val)
  {
    df_t_Ti$diff_t <- df_t_Ti %>% group_by(pid) %>% mutate(diff_t_Ti = exp(beta1*alpha1*t_Ti) - exp(beta1*alpha1*lag(t_Ti))) %>% pull(diff_t_Ti)
    
    df_t_Ti$Wi <- Wi
    
    df_t_Ti <- df_t_Ti[complete.cases(df_t_Ti),]
    
    l <- df_t_Ti %>% group_by(pid) %>% mutate(ll =  t_Ti^(gamma_h01-1)*exp(beta1*Wi)*diff_t/(beta1*alpha1)) %>% summarise(l = sum(ll), .groups = 'drop') %>% pull(l)
    
    return(-(sum(deltai*val*W[,1] - gamma_h01*gamma_h02^(-gamma_h01)*exp(val*W[,1] + beta1*(ui + alpha0))*l) + dnorm(val, mean = 0, sd = 100, log = T)))
  }
  
  #gam1_post <- stats4::mle(p_gam1_mle, start = list(val = gam1), method = "Brent", lower = -4, upper = 4)@coef
  
  #gam1_star <- rnorm(1, mean = gam1_post, sd = sd_gam1_star)
  gam1_star <- rnorm(1, mean = gam1, sd = sd_gam1_star)
  u_gam1 <- runif(1)
  
  log_r_gam1 <- p_gam1(gam1_star) + dnorm(gam1, mean = gam1_star, sd = sd_gam1_star, log = T) - p_gam1(gam1) - dnorm(gam1_star, mean = gam1, sd = sd_gam1_star, log = T)
  
  if(min(exp(log_r_gam1), 1) >= u_gam1)
  {
    update_gam1 <- 1
    return(list(gam1_star, update_gam1))
  } else {
    update_gam1 <- 0
    return(list(gam1, update_gam1))
  }
}


update_gam1 <- c()

gam1_i_1 <- gam1
gamma_h01_i_1 <- gamma_h01
gamma_h02_i_1 <- gamma_h02
beta1_i_1 <- beta1


system.time(for(i in 2:niter)
{
  res_gamma_h01 <- f_gamma_h01(ui = ui, alpha0 = alpha0, alpha1 = alpha1, gamma_h01 = gamma_h01_i_1, gamma_h02 = gamma_h02_i_1, beta1 = beta1_i_1, gam1 = gam1_i_1, Wi = Wi, df_t_Ti = df_t_Ti, y = df$FEV1, ni = ni, N = N, X = df[1:6], W = s_df[c(3:6, 2)], deltai = s_df$event, sd_gamma_h01_star = 0.1)
  
  gamma_h01_i <- res_gamma_h01[[1]]
  update_gamma_h01_i <- res_gamma_h01[[2]]
  
  res_gamma_h02 <- f_gamma_h02(ui = ui, alpha0 = alpha0, alpha1 = alpha1, gamma_h01 = gamma_h01_i, gamma_h02 = gamma_h02_i_1, beta1 = beta1_i_1, gam1 = gam1_i_1, Wi = Wi, df_t_Ti = df_t_Ti, y = df$FEV1, ni = ni, N, X = df[1:6], W = s_df[c(3:6, 2)], deltai = s_df$event, sd_gamma_h02_star = 0.08)
  
  gamma_h02_i <- res_gamma_h02[[1]]
  update_gamma_h02_i <- res_gamma_h02[[2]]
  
  res_beta1 <- f_beta1(ui = ui, alpha0 = alpha0, alpha1 = alpha1, gamma_h01 = gamma_h01_i, gamma_h02 = gamma_h02_i, beta1 = beta1_i_1, gam1 = gam1_i_1, Wi = Wi, df_t_Ti = df_t_Ti, y = df$FEV1, ni = ni, N = N, X = df[1:6], W = s_df[c(3:6, 2)], deltai = s_df$event, sd_beta1_star = 0.08)
  
  beta1_i <- res_beta1[[1]]
  update_beta1_i <- res_beta1[[2]]
  
  res_gam1 <- f_gam1(ui = ui, alpha0 = alpha0, alpha1 = alpha1, gamma_h01 = gamma_h01_i, gamma_h02 = gamma_h02_i, beta1 = beta1_i, gam1 = gam1_i_1,  Wi = Wi, df_t_Ti = df_t_Ti, y = df$FEV1, ni = ni, N = N, X = df[1:6], W = s_df[3:6], deltai = s_df$event, sd_gam1_star = 0.2)
  
  gam1_i <- res_gam1[[1]]
  update_gam1_i <- res_gam1[[2]]
  
  if(i == 2 | i %% 100 == 0){
    print(paste0("iteration: ", i, " (",i*100/niter,"%)"))
    
    print(paste0("gamma_h01[",i,"] = ", gamma_h01_i, "-- AR:", round(prop.table(table(update_gamma_h01))[2], 2)))
    print(paste0("gamma_h02[",i,"] = ", gamma_h02_i, "-- AR:", round(prop.table(table(update_gamma_h02))[2], 2)))
    print(paste0("beta1[",i,"] = ", beta1_i, "-- AR:", round(prop.table(table(update_beta1))[2], 2)))
    print(paste0("gam1[",i,"] = ", gam1_i, "-- AR:", round(prop.table(table(update_gam1))[2], 2)))
    
    print("------------------------------------------------")
  }
  
  
  if(i %% lag == 0)
  {
    gam1[i/lag] <- gam1_i
    
    gamma_h01[i/lag] <- gamma_h01_i
    gamma_h02[i/lag] <- gamma_h02_i
    beta1[i/lag] <- beta1_i
    
    update_gam1[i/lag] <- res_gam1[[2]]
    
    update_gamma_h01[i/lag] <- res_gamma_h01[[2]]
    update_gamma_h02[i/lag] <- res_gamma_h02[[2]]
    
    update_beta1[i/lag] <- res_beta1[[2]]
  }
  gam1_i_1 <- gam1_i
  
  gamma_h01_i_1 <- gamma_h01_i
  gamma_h02_i_1 <- gamma_h02_i
  beta1_i_1 <- beta1_i
})

parameters <- data.frame(gamma_h01, gamma_h02, beta1, gam1)

saveRDS(parameters, paste0("survivalSM_simulation", as.numeric(args[1]), "_chain2.rds"))

updates <- data.frame(update_gamma_h01,update_gamma_h02, update_beta1, update_gam1)

saveRDS(updates, paste0("survivalSM_simulation", as.numeric(args[1]), "_updates_chain2.rds"))
