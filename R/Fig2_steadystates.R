# This script creates the surface plots used in Figure 2.
library(foreach)
library(parallel)
library(doParallel)
library(reshape2)
library(dplyr)

source("R/def_pars.R")
source("R/init_env.R")
source("R/run_coral.R")


# Set time vector and default parameters
time <- seq(1,120,1)
pars <- def_pars()
pars$jCO2a <- 1  # having higher carbon delivery from host reduces S:H ratios and causes photosynthesis to be light-limited more frequently.

# Set values of L and N at which to get steady state values
at <- expand.grid(L=seq(from=0, to=40, length.out=41),
                  N=seq(from=0, to=2e-6, length.out=41))

# Set up cluster for parallel processing
cl <- makeCluster(detectCores())  # Initiate cluster
registerDoParallel(cl)

ptm <- proc.time()
# Run simulations in parallel for each combination of L and N values
steady_states <- foreach(i=1:nrow(at), .combine='rbind', .packages='dplyr') %dopar% {
  L <- at[i,]$L; N <- at[i,]$N
  env <- init_env(time=time, L=c(L,L,0), N=c(N,N,0), X=c(0,0,0))
  run <- run_coral(time=time, env=env, pars=pars)
  ss <- lapply(run[c("H", "S")], function(x) x[length(time), ])
  gr <- ss$H$dH.Hdt
  sh <- ss$S$S/ss$H$H
  hl <- with(ss, log(  pmin(S$rhoC*S$S/H$H + H$jX, pars$jHGm) / pmin((H$jN + pars$nNX*H$jX + H$rNH) / pars$nNH, pars$jHGm)  ))
  sl <- with(ss, log(  pmin(S$jCP, pars$jSGm)   /   pmin((H$rhoN*H$H/S$S + S$rNS)/pars$nNS, pars$jSGm)  ))
  ee <- with(ss, max(0, S$jL - (S$jCP/run$pars$nLC + run$pars$jNPQ)))
  pl <- with(ss, log(   pmin((H$jCO2 + H$rCH)*H$H/S$S + S$rCS, pars$jCPm)   /   pmin(S$jL * pars$nLC, pars$jCPm)    ))
  data.frame(gr=gr, sh=sh, hl=hl, sl=sl, ee=ee, pl=pl)
}
ss <- cbind(at, steady_states)
save(ss, file = "output/ss.Rdata")

stopCluster(cl)  # Stop cluster
ptm - proc.time()  # 2601 simulations with 100 time steps each = 427 sec
                   # 10201 simulations with 100 time steps each = 1670 sec






