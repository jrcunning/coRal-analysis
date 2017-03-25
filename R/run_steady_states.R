run_steady_states <- function(pars, at, outfile=NULL, runtime=F, food=0) {
  
  # Load libraries
  require(foreach)
  require(parallel)
  require(doParallel)
  require(reshape2)
  require(dplyr)
  require(coRal)
  
  # Set up cluster for parallel processing
  cl <- makeCluster(detectCores())  # Initiate cluster
  registerDoParallel(cl)
  
  # Run simulations in parallel for each combination of "at" values
  start <- proc.time()  # Start timer...
  # Calculate steady states...
  steady_states <- foreach(i=1:nrow(at), .combine='rbind', .packages='dplyr') %dopar% {
    env <- list(L=at[i,]$L, N=at[i,]$N, X=food)
    run <- coRal::run_coral_ss(env=env, pars=pars, dt=0.1)
    ss <- lapply(run[c("H", "S")], function(x) x[nrow(x), ])
    gr <- ss$H$dH.Hdt
    sh <- ss$S$S/ss$H$H
    data.frame(gr=gr, sh=sh)
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
