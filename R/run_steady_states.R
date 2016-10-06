run_steady_states <- function(time, pars, at, outfile=NULL, runtime=F) {
  
  # Load libraries
  require(foreach)
  require(parallel)
  require(doParallel)
  require(reshape2)
  require(dplyr)
  
  # Load functions
  source("R/init_env.R", local=TRUE)
  source("R/run_coral.R", local=TRUE)
  
  # Set up cluster for parallel processing
  cl <- makeCluster(detectCores())  # Initiate cluster
  registerDoParallel(cl)
  
  # Run simulations in parallel for each combination of "at" values
  start <- proc.time()  # Start timer...
  # Calculate steady states...
  steady_states <- foreach(i=1:nrow(at), .combine='rbind', .packages='dplyr') %dopar% {
    L <- at[i,]$L; N <- at[i,]$N
    env <- init_env(time=time, L=c(L,L,0), N=c(N,N,0), X=c(0,0,0))
    run <- run_coral(time=time, env=env, pars=pars)
    ss <- lapply(run[c("H", "S")], function(x) x[length(time), ])
    gr <- ss$H$dH.Hdt
    sh <- ss$S$S/ss$H$H
    hl <- with(ss, log(  pmin(S$rhoC*S$S/H$H + H$jX, pars$jHGm) / pmin((H$jN + pars$nNX*H$jX + H$rNH) / pars$nNH, pars$jHGm)  ))
    sl <- with(ss, log(  pmin(S$jCP, pars$jSGm)   /   pmin((H$rhoN*H$H/S$S + S$rNS)/pars$nNS, pars$jSGm)  ))
    ee <- with(ss, max(0, S$jL - (S$jCP/run$pars$nLC + run$pars$kNPQ)))
    pl <- with(ss, log(   pmin((H$jCO2 + H$rCH)*H$H/S$S + S$rCS, pars$jCPm)   /   pmin(S$jL * pars$nLC, pars$jCPm)    ))
    data.frame(gr=gr, sh=sh, hl=hl, sl=sl, ee=ee, pl=pl)
  }
  stopCluster(cl)  # Stop cluster
  if (runtime==T) print(proc.time() - start)  # Print time elapsed if runtime==T
  
  # Collect results in ss
  ss <- cbind(at, steady_states)
  
  # If output filename give, save ss as Rdata
  if (!is.null(outfile)) save(ss, file=outfile)
  
  # Return steady state values
  return(ss)
}
