#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(invgamma)
library(dplyr)
library(mvtnorm)
library(stats4)
library(pracma)

set.seed(12345 + as.numeric(args[1]))

df1 <- read.csv(paste0("simudata", as.numeric(args[1]), ".csv"))

niter <- 100000

N = nrow(df1)
Ni = length(unique(df1$pid))

y = df1$FEV1

ni <- aggregate(df1$lamdxage, by = list(df1$pid), length)$x

df <- df1[c(2, 3, 8, 5, 9, 10, 11)]

ni_l = c(1, (cumsum(ni)+1)[-length(ni)])
ni_u = cumsum(ni)

pt = unique(df$pid)

alpha0 =  c()
alpha1 = c()
sigmasq = c()
tausq = c()
H = c()
kappa = c()
ui <- matrix(NA, ncol = niter/25, nrow = Ni)
Wi <- matrix(NA, ncol = niter/25, nrow = N)

alpha0[1] =  as.numeric(args[2])
alpha1[1] = as.numeric(args[3])
sigmasq[1] = as.numeric(args[4])
tausq[1] = as.numeric(args[5])
H[1] <- as.numeric(args[6])
kappa[1] <- as.numeric(args[7])
ui[,1] = readRDS(paste0("ui", as.numeric(args[1]),".rds"))
Wi[,1] <- readRDS(paste0("Wi", as.numeric(args[1]),".rds"))

f_sigmasq <- function(ui, alpha0, alpha1, Wi, y, ni, N, X)
{
  rinvgamma(1, shape = N/2 + 0.001, rate = sum((y - rep(ui, ni) - alpha0 - alpha1*X[,2] - Wi)^2)/2 + 0.001)
}

f_tausq <- function(ui, Ni)
{
  rinvgamma(1, shape = Ni/2 + 0.001, rate = sum(ui^2)/2 + 0.001)
}

f_kappa <- function(Wij, H, N, Ni, ni, ni_l, ni_u, t)
{
  f_kappa_rate <- 0
  for(i in 1:Ni)
  {
    Ci <- matrix(NA, ni[i], ni[i])
    t_Ti <- t[ni_l[i]:ni_u[i]]
    Wi <- Wij[ni_l[i]:ni_u[i]]
    
    for(j in 1:(ni[i]-1)){
      Ci[j, j] = t_Ti[j]^(2*H + 2)/(2*H + 2)
      for(k in (j+1):ni[i]){
        Ci[j, k] = t_Ti[j]^(2*H + 2)/(2*H + 2) + 0.5*(t_Ti[j]*t_Ti[k]^(2*H + 1) + t_Ti[k]*t_Ti[j]^(2*H + 1) + (t_Ti[k] - t_Ti[j])^(2*H + 2)/(2*H + 2) - 2*t_Ti[j]^(2*H + 2) - (t_Ti[k]^(2*H + 2) - t_Ti[j]^(2*H + 2))/(2*H + 2))/(2*H + 1)
        Ci[k, j] = Ci[j, k]
      }
    } 
    Ci[ni[i], ni[i]] = t_Ti[ni[i]]^(2*H + 2)/(2*H + 2)
    
    if(isposdef(Ci))
    {
      Ci = Ci
    } else {
      Ci = Ci + 1e-11*diag(ni[i])
    }
    
    f_kappa_rate <- f_kappa_rate + 0.5*t(Wi)%*%solve(Ci)%*%t(t(Wi))
  }
  rinvgamma(1, shape = N/2 + 0.001, rate = f_kappa_rate + 0.001)
}

f_H <- function(kappa, H, Wij, ni, ni_l, ni_u, Ni, t, sd_H_star)
{
  H_current <- log(H/(1-H))
  p_H_mle <- function(val)
  {
    H_llike <- 0
    for(i in 1:Ni)
    {
      t_Ti <- t[ni_l[i]:ni_u[i]]
      Ci <- matrix(NA, ni[i], ni[i])
      
      for(j in 1:(ni[i]-1)){
        Ci[j, j] = t_Ti[j]^(2*(exp(val)/(1 + exp(val))) + 2)/(2*(exp(val)/(1 + exp(val))) + 2)
        for(k in (j+1):ni[i]){
          Ci[j, k] = t_Ti[j]^(2*(exp(val)/(1 + exp(val))) + 2)/(2*(exp(val)/(1 + exp(val))) + 2) + 0.5*(t_Ti[j]*t_Ti[k]^(2*(exp(val)/(1 + exp(val))) + 1) + t_Ti[k]*t_Ti[j]^(2*(exp(val)/(1 + exp(val))) + 1) + (t_Ti[k] - t_Ti[j])^(2*(exp(val)/(1 + exp(val))) + 2)/(2*(exp(val)/(1 + exp(val))) + 2) - 2*t_Ti[j]^(2*(exp(val)/(1 + exp(val))) + 2) - (t_Ti[k]^(2*(exp(val)/(1 + exp(val))) + 2) - t_Ti[j]^(2*(exp(val)/(1 + exp(val))) + 2))/(2*(exp(val)/(1 + exp(val))) + 2))/(2*(exp(val)/(1 + exp(val))) + 1)
          Ci[k, j] = Ci[j, k]
        }
      } 
      Ci[ni[i], ni[i]] = t_Ti[ni[i]]^(2*(exp(val)/(1 + exp(val))) + 2)/(2*(exp(val)/(1 + exp(val))) + 2)
      
      Wi <- Wij[ni_l[i]:ni_u[i]]
      
      if(isposdef(kappa*Ci))
      {
        Ci = Ci
      } else {
        Ci = Ci + 1e-11*diag(ni[i])
      }
      
      H_llike <- H_llike - 0.5*(log(det(Ci)) + t(Wi)%*%solve(kappa*Ci)%*%t(t(Wi)))
    }
    
    return(-(H_llike + val - 2*log(1 + exp(val)) + dnorm(val, mean = 0, sd = 100, log = T)))
  }
  
  p_H_hessian <- function(val)
  {
    H_llike <- 0
    for(i in 1:Ni)
    {
      t_Ti <- t[ni_l[i]:ni_u[i]]
      Ci <- matrix(NA, ni[i], ni[i])
      
      for(j in 1:(ni[i]-1)){
        Ci[j, j] = t_Ti[j]^(2*(exp(val)/(1 + exp(val))) + 2)/(2*(exp(val)/(1 + exp(val))) + 2)
        for(k in (j+1):ni[i]){
          Ci[j, k] = t_Ti[j]^(2*(exp(val)/(1 + exp(val))) + 2)/(2*(exp(val)/(1 + exp(val))) + 2) + 0.5*(t_Ti[j]*t_Ti[k]^(2*(exp(val)/(1 + exp(val))) + 1) + t_Ti[k]*t_Ti[j]^(2*(exp(val)/(1 + exp(val))) + 1) + (t_Ti[k] - t_Ti[j])^(2*(exp(val)/(1 + exp(val))) + 2)/(2*(exp(val)/(1 + exp(val))) + 2) - 2*t_Ti[j]^(2*(exp(val)/(1 + exp(val))) + 2) - (t_Ti[k]^(2*(exp(val)/(1 + exp(val))) + 2) - t_Ti[j]^(2*(exp(val)/(1 + exp(val))) + 2))/(2*(exp(val)/(1 + exp(val))) + 2))/(2*(exp(val)/(1 + exp(val))) + 1)
          Ci[k, j] = Ci[j, k]
        }
      } 
      Ci[ni[i], ni[i]] = t_Ti[ni[i]]^(2*(exp(val)/(1 + exp(val))) + 2)/(2*(exp(val)/(1 + exp(val))) + 2)
      
      Wi <- Wij[ni_l[i]:ni_u[i]]
      
      if(isposdef(kappa*Ci))
      {
        Ci = Ci
      } else {
        Ci = Ci + 1e-11*diag(ni[i])
      }
      
      H_llike <- H_llike - 0.5*(log(det(Ci)) + t(Wi)%*%solve(kappa*Ci)%*%t(t(Wi)))
    }
    
    return(H_llike + val - 2*log(1 + exp(val)) + dnorm(val, mean = 0, sd = 100, log = T))
  }
  
  #H1 <- stats4::mle(p_H_mle, start = list(val = H_current))@coef
  #H1_star <- rnorm(1, mean = H1, sd = sd_H_star)
  
  #log_r_H <- p_H_hessian(H1_star) + dnorm(H_current, mean = H1, sd = sd_H_star, log = T) - p_H_hessian(H_current) - dnorm(H1_star, mean = H1, sd = sd_H_star, log = T)
  
  H1_star <- rnorm(1, mean = H_current, sd = sd_H_star)
  log_r_H <- p_H_hessian(H1_star) + dnorm(H_current, mean = H1_star, sd = sd_H_star, log = T) - p_H_hessian(H_current) - dnorm(H1_star, mean = H_current, sd = sd_H_star, log = T)
  
  r_H <- exp(log_r_H)
  
  u_H <- runif(1)
  
  if(min(r_H, 1) >= u_H){
    
    update_H <- 1
    return(list(exp(H1_star)/(1 + exp(H1_star)), update_H))
  } else 
  {
    update_H <- 0
    return(list(exp(H_current)/(1 + exp(H_current)), update_H))
  }
}

update_H <- c()

f_ui <- function(alpha0, alpha1, tausq, sigmasq, Wij, yij, Xij, ni)
{
  return(rnorm(1, mean =  tausq*sum(yij - alpha0 - alpha1*Xij[,2] - Wij)/(tausq*ni + sigmasq), sd = sqrt((ni/sigmasq + 1/tausq)^(-1))))
}

f_Wi <- function(ui, alpha0, alpha1, sigmasq, H, kappa, y, X, ni)
{
  Ci <- matrix(NA, ni, ni)
  t_Ti <- X$timeSince_t0
  
  for(j in 1:(ni-1)){
    Ci[j, j] = t_Ti[j]^(2*H + 2)/(2*H + 2)
    for(k in (j+1):ni){
      Ci[j, k] = t_Ti[j]^(2*H + 2)/(2*H + 2) + 0.5*(t_Ti[j]*t_Ti[k]^(2*H + 1) + t_Ti[k]*t_Ti[j]^(2*H + 1) + (t_Ti[k] - t_Ti[j])^(2*H + 2)/(2*H + 2) - 2*t_Ti[j]^(2*H + 2) - (t_Ti[k]^(2*H + 2) - t_Ti[j]^(2*H + 2))/(2*H + 2))/(2*H + 1)
      Ci[k, j] = Ci[j, k]
    }
  } 
  Ci[ni, ni] = t_Ti[ni]^(2*H + 2)/(2*H + 2)
  
  if(isposdef(kappa*Ci))
  {
    Ci = Ci
  } else {
    Ci = Ci + 1e-11*diag(ni)
  }
  
  return(mvtnorm::rmvnorm(1, mean = solve(diag(ni)/sigmasq + solve(kappa*Ci))%*%(diag(ni)%*%t(t(y - rep(ui, ni) - alpha0 - alpha1*X[,2])))/sigmasq, sigma = solve(diag(ni)/sigmasq + solve(kappa*Ci))))
}

f_alpha0 <- function(ui, alpha1, sigmasq, Wi, y, ni, N, X)
{
  rnorm(1, mean = 10^4*sum(y - rep(ui, ni) - alpha1*X[,2] - Wi)/(sigmasq + 10^4*N), sd = sqrt(10^4*sigmasq/(sigmasq + 10^4*N)))
}

update_alpha0 <- c()

f_alpha1 <- function(ui, alpha0, sigmasq, Wi, y, ni, N, X)
{
  return(rnorm(1, mean = 10^4*sum(X[,2]*(y - rep(ui, ni) - alpha0 - Wi))/(sigmasq + 10^4*sum(X[,2]^2)), sd = sqrt(10^4*sigmasq/(sigmasq + 10^4*sum(X[,2]^2)))))
}

update_alpha1 <- c()

Wi_i_1 <- Wi[,1]
ui_i_1 <- ui[,1]

ui_i <- c()
Wi_i <- c()

system.time(for(i in 2:niter)
{
  sigmasq[i] <- f_sigmasq(ui = ui_i_1, alpha0 = alpha0[i-1], alpha1 = alpha1[i-1], Wi = Wi_i_1, y = df$FEV1, ni = ni, N = N, X = df[1:6])
  
  tausq[i] <- f_tausq(ui = ui_i_1, Ni = Ni)
  
  kappa[i] <- f_kappa(Wij = Wi_i_1, H = H[i-1], N = N, Ni = Ni, ni = ni, ni_l = ni_l, ni_u = ni_u, t = df$timeSince_t0)
  
  res_H <- f_H(H = H[i - 1], kappa = kappa[i], Wij = Wi_i_1, ni = ni, ni_l = ni_l, ni_u = ni_u, Ni = Ni, t = df$timeSince_t0, sd_H_star = 0.2)
  
  H[i] <- res_H[[1]]
  update_H[i] <- res_H[[2]]
  
  for(j in 1:Ni)
  {
    df_i <- subset(df, pid == pt[j])
    ui_i[j] <- f_ui(alpha0 = alpha0[i-1], alpha1 = alpha1[i-1], tausq = tausq[i], sigmasq = sigmasq[i], Wij = Wi_i_1[ni_l[j]:ni_u[j]], yij = df_i$FEV1, Xij = df_i[1:6], ni = ni[j])
  }
  
  for(j in 1:Ni)
  {
    df_i <- subset(df, pid == pt[j])
    
    Wi_i[ni_l[j]:ni_u[j]] <- f_Wi(ui = ui_i[j], alpha0 = alpha0[i-1], alpha1= alpha1[i-1], sigmasq = sigmasq[i], kappa = kappa[i], H = H[i], y = df_i$FEV1, X = df_i[1:6], ni = ni[j])
  }
  
  alpha0[i] <- f_alpha0(ui = ui_i, alpha1 = alpha1[i-1], sigmasq = sigmasq[i], Wi = Wi_i, y = df$FEV1, ni = ni, N = N, X = df[1:6])
  
  alpha1[i] <- f_alpha1(ui = ui_i, alpha0 = alpha0[i], sigmasq = sigmasq[i], Wi = Wi_i, y = df$FEV1, ni = ni, N = N, X = df[1:6])
  
  if(i %% 25 == 0)
  {
    Wi[,i/25] <- Wi_i
    ui[,i/25] <- ui_i
  }
  
  if(i %% 1000 == 0)
  {
    saveRDS(data.frame(alpha0, alpha1, sigmasq, tausq, kappa, H), paste0("longitudinalSM_simulation", as.numeric(args[1]), "_chain1.rds"))
  }
  if(i %% 100 == 0)
  {
    print(paste0("iteration: ", i, "(",i*100/niter,"%)"))
    print(paste0("alpha0[",i,"] = ", alpha0[i]))
    print(paste0("alpha1[",i,"] = ", alpha1[i]))
    print(paste0("Kappa[",i,"] = ", kappa[i]))
    print(paste0("H[",i,"] = ", H[i], "-- Acceptance Prob: ", round(prop.table(table(update_H))[2], 2)))
    print(paste0("sigmasq[",i,"] = ", sigmasq[i]))
    print(paste0("tausq[",i,"] = ", tausq[i]))
    
    print("------------------------------------------------")
  }
  
  Wi_i_1 <- Wi_i
  ui_i_1 <- ui_i
})

parameters <- data.frame(alpha0, alpha1, sigmasq, tausq, kappa, H)
updates <- update_H
random_terms <- data.frame(t(Wi), t(ui))

saveRDS(parameters, paste0("longitudinalSM_simulation", as.numeric(args[1]), "_chain1.rds"))
saveRDS(random_terms, paste0("longitudinalSM_simulation", as.numeric(args[1]), "_randomterms_chain1.rds"))
saveRDS(updates, paste0("longitudinalSM_simulation", as.numeric(args[1]), "_updates_chain1.rds"))
