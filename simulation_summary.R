library(coda)

TrueValues <- readRDS("TrueValues.rds")

library(coda)

### Part 1###

for(k in 1:100)
{
  df1 <- read.csv(paste0("simudata", as.numeric(k), ".csv"))
  
  N = nrow(df1)
  Ni = length(unique(df1$pid))
  ni <- aggregate(df1$lamdxage, by = list(df1$pid), length)$x
  
  chain1 <- readRDS(paste0("longitudinalSM_simulation",k,"_randomterms_chain1.rds"))
  chain2 <- readRDS(paste0("longitudinalSM_simulation",k,"_randomterms_chain2.rds"))
  
  summary_re_stoc <- matrix(NA, ncol = 13, nrow = N + Ni)
  for(i in 1:dim(chain1)[2])
  {
    para_mcmc <- mcmc.list(as.mcmc(chain1[3600:4000, i]), as.mcmc(chain2[3600:4000, i]))
    
    summary_re_stoc[i,] <- c(drop(summary(para_mcmc)[[1]]), drop(summary(para_mcmc)[[2]]), drop(gelman.diag(para_mcmc)[[1]]), c(geweke.diag(para_mcmc)[[1]])$z, c(geweke.diag(para_mcmc)[[2]])$z)
  }
  
  summary_re_stoc <- as.data.frame(summary_re_stoc)
  
  colnames(summary_re_stoc) <- c("Mean", "SD", "Naive SE", "Time-series SE", "2.5% quantile", "25% quantile", "50% quantile", "75% quantile", "97.5% quantile", "BGR diagnostics", "Upper C.I. of BGR diagnostics", "Geweke's convergence diagnostics (chain 1)", "Geweke's convergence diagnostic (chain 2)")
  
  j <- c()
  for(i in 1:Ni)
  {
    j <- c(j, seq(1:ni[i]))
  }
  
  rownames(summary_re_stoc) <- c(paste0("W", rep(1:Ni, times = ni),"," ,j), paste0("U", seq(1:Ni)))
  
  saveRDS(summary_re_stoc, paste0("longitudinalSM_random_terms_estimates_simulation",k,".rds"))

  main_chain1 <- readRDS(paste0("longitudinalSM_simulation", k,"_chain1.rds"))
  main_chain2 <- readRDS(paste0("longitudinalSM_simulation", k,"_chain2.rds"))
  
  main_para_diagnostics <- matrix(NA, ncol = 13, nrow = 6)
  
  alpha0_mcmc <- mcmc.list(as.mcmc(main_chain1[50000:100000, "alpha0"]), as.mcmc(main_chain2[50000:100000, "alpha0"]))
  
  main_para_diagnostics[1,] <- c(drop(gelman.diag(alpha0_mcmc)[[1]]), c(geweke.diag(alpha0_mcmc)[[1]])$z, c(geweke.diag(alpha0_mcmc)[[2]])$z, drop(summary(alpha0_mcmc)[[1]]), drop(summary(alpha0_mcmc)[[2]]))
  
  alpha1_mcmc <- mcmc.list(as.mcmc(main_chain1[50000:100000, "alpha1"]), as.mcmc(main_chain2[50000:100000, "alpha1"]))
  
  main_para_diagnostics[2,] <- c(drop(gelman.diag(alpha1_mcmc)[[1]]), c(geweke.diag(alpha1_mcmc)[[1]])$z, c(geweke.diag(alpha1_mcmc)[[2]])$z, drop(summary(alpha1_mcmc)[[1]]), drop(summary(alpha1_mcmc)[[2]]))
  
  sigmasq_mcmc <- mcmc.list(as.mcmc(main_chain1[50000:100000, "sigmasq"]), as.mcmc(main_chain2[50000:100000, "sigmasq"]))
  
  main_para_diagnostics[3,] <- c(drop(gelman.diag(sigmasq_mcmc)[[1]]), c(geweke.diag(sigmasq_mcmc)[[1]])$z, c(geweke.diag(sigmasq_mcmc)[[2]])$z, drop((summary(sigmasq_mcmc)[[1]])), drop(summary(sigmasq_mcmc)[[2]]))
  
  tausq_mcmc <- mcmc.list(as.mcmc(main_chain1[50000:100000, "tausq"]), as.mcmc(main_chain2[50000:100000, "tausq"]))
  
  main_para_diagnostics[4,] <- c(drop(gelman.diag(tausq_mcmc)[[1]]), c(geweke.diag(tausq_mcmc)[[1]])$z, c(geweke.diag(tausq_mcmc)[[2]])$z, drop(summary(tausq_mcmc)[[1]]), drop(summary(tausq_mcmc)[[2]]))
  
  kappa_mcmc <- mcmc.list(as.mcmc(main_chain1[50000:100000, "kappa"]), as.mcmc(main_chain2[50000:100000, "kappa"]))
  
  main_para_diagnostics[5,] <- c(drop(gelman.diag(kappa_mcmc)[[1]]), c(geweke.diag(kappa_mcmc)[[1]])$z, c(geweke.diag(kappa_mcmc)[[2]])$z, drop(summary(kappa_mcmc)[[1]]), drop(summary(kappa_mcmc)[[2]]))
  
  H_mcmc <- mcmc.list(as.mcmc(main_chain1[50000:100000, "H"]), as.mcmc(main_chain2[50000:100000, "H"]))
  
  main_para_diagnostics[6,] <- c(drop(gelman.diag(H_mcmc)[[1]]), c(geweke.diag(H_mcmc)[[1]])$z, c(geweke.diag(H_mcmc)[[2]])$z, drop(summary(H_mcmc)[[1]]), drop(summary(H_mcmc)[[2]]))
  
  main_para_diagnostics <- as.data.frame(main_para_diagnostics)
  
  colnames(main_para_diagnostics) <- c("BGR diagnostics", "Upper C.I. of BGR diagnostics", "Geweke's convergence diagnostics (chain 1)", "Geweke's convergence diagnostic (chain 2)", "Mean", "SD", "Naive SE", "Time-series SE", "2.5% quantile", "25% quantile", "50% quantile", "75% quantile", "97.5% quantile")
  
  rownames(main_para_diagnostics) <- c("alpha0", "alpha1", "sigmasq", "tausq", "kappa", "H")
  
  saveRDS(main_para_diagnostics, paste0("longitudinalSM_mainpara_estimates_simulation",k,".rds"))
 
  png(paste0("Traceplots_longitudinalSM_Simulation", k, ".png"), height = 1200, width = 1200, pointsize = 20)
  par(mfrow = c(3, 2))
  plot(main_chain1[,"alpha0"], type = "l", ylab = "alpha0")
  lines(main_chain2[,"alpha0"], col = "light green")
  abline(h = TrueValues[1, 1], col = "red", lty = 2, lwd = 2)
  
  plot(main_chain1[,"alpha1"], type = "l", ylab = "alpha1")
  lines(main_chain2[,"alpha1"], col = "light green")
  abline(h = TrueValues[1, 2], col = "red", lty = 2, lwd = 2)
  
  plot(main_chain1[,"sigmasq"], type = "l", ylab = "sigmasq")
  lines(main_chain2[,"sigmasq"], col = "light green")
  abline(h = TrueValues[1, 3]^2, col = "red", lty = 2, lwd = 2)
  
  plot(main_chain1[,"tausq"], type = "l", ylab = "tausq")
  lines(main_chain2[,"tausq"], col = "light green")
  abline(h = TrueValues[1, 4]^2, col = "red", lty = 2, lwd = 2)
  
  plot(main_chain1[,"kappa"], type = "l", ylab = "kappa")
  lines(main_chain2[,"kappa"], col = "light green")
  abline(h = TrueValues[1, 5], col = "red", lty = 2, lwd = 2)
  
  plot(main_chain1[,"H"], type = "l", ylab = "H")
  lines(main_chain2[,"H"], col = "light green")
  abline(h = TrueValues[1, 6], col = "red", lty = 2, lwd = 2)
  mtext(paste0("Simulation", i),side=3,outer=TRUE,padj=3)
  dev.off() 
}

### Part 2###

for(i in 1:100)
{
  chain1 <- readRDS(paste0("survivalSM_simulation",i,"_chain1.rds"))
  chain2 <- readRDS(paste0("survivalSM_simulation",i,"_chain2.rds"))

  png(paste0("Traceplots_survivalSM_Simulation", i, ".png"), height = 800, width = 1200, pointsize = 20)
  par(mfrow = c(2, 2))
  plot(chain1[,"gamma_h01"], type = "l", ylab = "gamma_h01")
  lines(chain2[,"gamma_h01"], col = "light green")
  abline(h = TrueValues[1, 9], col = "red", lty = 2, lwd = 2)
  
  plot(chain1[,"gamma_h02"], type = "l", ylab = "gamma_h02")
  lines(chain2[,"gamma_h02"], col = "light green")
  abline(h = TrueValues[1, 10], col = "red", lty = 2, lwd = 2)
  
  plot(chain1[,"gam1"], type = "l", ylab = "gam1")
  lines(chain2[,"gam1"], col = "light green")
  abline(h = TrueValues[1, 7], col = "red", lty = 2, lwd = 2)
  
  plot(chain1[,"beta1"], type = "l", ylab = "beta1")
  lines(chain2[,"beta1"], col = "light green")
  abline(h = TrueValues[1, 8], col = "red", lty = 2, lwd = 2)
  dev.off()
  
  main_para_diagnostics <- matrix(NA, ncol = 11, nrow = 4)
  
  gamma_h01_mcmc <- mcmc.list(as.mcmc(chain1[25000:50000, "gamma_h01"]), as.mcmc(chain2[25000:50000, "gamma_h01"]))
  
  main_para_diagnostics[1,] <- c(TrueValues[1, 9], c(drop(summary(gamma_h01_mcmc, quantiles = c(0.025, 0.975))[[1]]), drop(summary(gamma_h01_mcmc, quantiles = c(0.025, 0.975))[[2]]), drop(gelman.diag(gamma_h01_mcmc)[[1]]), c(geweke.diag(gamma_h01_mcmc)[[1]])$z, c(geweke.diag(gamma_h01_mcmc)[[2]])$z))
  
  gamma_h02_mcmc <- mcmc.list(as.mcmc(chain1[25000:50000, "gamma_h02"]), as.mcmc(chain2[25000:50000, "gamma_h02"]))
  
  main_para_diagnostics[2,] <- c(TrueValues[1, 10], c(drop(summary(gamma_h02_mcmc, quantiles = c(0.025, 0.975))[[1]]), drop(summary(gamma_h02_mcmc, quantiles = c(0.025, 0.975))[[2]]), drop(gelman.diag(gamma_h02_mcmc)[[1]]), c(geweke.diag(gamma_h02_mcmc)[[1]])$z, c(geweke.diag(gamma_h02_mcmc)[[2]])$z))
  
  beta1_mcmc <- mcmc.list(as.mcmc(chain1[25000:50000, "beta1"]), as.mcmc(chain2[25000:50000, "beta1"]))
  
  main_para_diagnostics[3,] <- c(TrueValues[1, 8], c(drop(summary(beta1_mcmc, quantiles = c(0.025, 0.975))[[1]]), drop(summary(beta1_mcmc, quantiles = c(0.025, 0.975))[[2]]), drop(gelman.diag(beta1_mcmc)[[1]]), c(geweke.diag(beta1_mcmc)[[1]])$z, c(geweke.diag(beta1_mcmc)[[2]])$z))
  
  gam1_mcmc <- mcmc.list(as.mcmc(chain1[25000:50000, "gam1"]), as.mcmc(chain2[25000:50000, "gam1"]))
  
  main_para_diagnostics[4,] <- c(TrueValues[1, 7], c(drop(summary(gam1_mcmc, quantiles = c(0.025, 0.975))[[1]]), drop(summary(gam1_mcmc, quantiles = c(0.025, 0.975))[[2]]), drop(gelman.diag(gam1_mcmc)[[1]]), c(geweke.diag(gam1_mcmc)[[1]])$z, c(geweke.diag(gam1_mcmc)[[2]])$z))
  
  main_para_diagnostics <- as.data.frame(main_para_diagnostics)
  
  colnames(main_para_diagnostics) <- c("True value","Mean", "SD", "Naive SE", "Time-series SE", "2.5% quantile", "97.5% quantile", "BGR diagnostics", "Upper C.I. of BGR diagnostics", "Geweke's convergence diagnostics (chain 1)", "Geweke's convergence diagnostic (chain 2)")
  
  rownames(main_para_diagnostics) <- c("gamma_h01", "gamma_h02", "beta1", "gamma1")
  
  saveRDS(main_para_diagnostics, paste0("main_para_estimates_simulation", i, ".rds"))
  
}

