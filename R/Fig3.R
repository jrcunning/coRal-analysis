# Load functions
sapply(c("R/def_pars.R",
         "R/run_coral_ss.R", 
         "R/sensitivity.R"), 
       source, .GlobalEnv)

# Set up time vector and parameters
defpars <- def_pars()  # Get default parameters
defpars <- replace(defpars, "initS", 1)  # Start with high S biomass to avoid alternate steady state of negative growth

# Define environments for conducting sensitivity analysis
envs <- list(
  #LLLN=list(L=2, N=1e-7, X=0),
  #LLHN=list(L=2, N=4e-6, X=0),
  HLLN=list(L=15, N=1e-7, X=0)#,
  #HLHN=list(L=20, N=4e-6, X=0)
)

# Define relative changes in each parameter to evaluate
levs <- seq(0.5, 1.5, 0.05)  # Set levels of % change in parameter value to evaluate

# Define all simulations to be run for a given parameter
sims <- expand.grid(env=names(envs), change=levs)
n <- nrow(sims)

# Set up cluster for parallel processing
cl <- makeCluster(detectCores())  # Initiate cluster
registerDoParallel(cl)

# Run sensitivity analysis in each environment for each parameter
jHT0.sens <- foreach(i=1:n, .combine=rbind) %dopar% with(sims[i,], sens(env=envs[[env]], pars=defpars, par="jHT0", change=change))
jNm.sens <- foreach(i=1:n, .combine=rbind) %dopar% with(sims[i,], sens(env=envs[[env]], pars=defpars, par="jNm", change=change))
jHGm.sens <- foreach(i=1:n, .combine=rbind) %dopar% with(sims[i,], sens(env=envs[[env]], pars=defpars, par="jHGm", change=change))
kCO2.sens <- foreach(i=1:n, .combine=rbind) %dopar% with(sims[i,], sens(env=envs[[env]], pars=defpars, par="kCO2", change=change))
KN.sens <- foreach(i=1:n, .combine=rbind) %dopar% with(sims[i,], sens(env=envs[[env]], pars=defpars, par="KN", change=change))
jST0.sens <- foreach(i=1:n, .combine=rbind) %dopar% with(sims[i,], sens(env=envs[[env]], pars=defpars, par="jST0", change=change))
kNPQ.sens <- foreach(i=1:n, .combine=rbind) %dopar% with(sims[i,], sens(env=envs[[env]], pars=defpars, par="kNPQ", change=change))
kROS.sens <- foreach(i=1:n, .combine=rbind) %dopar% with(sims[i,], sens(env=envs[[env]], pars=defpars, par="kROS", change=change))
astar.sens <- foreach(i=1:n, .combine=rbind) %dopar% with(sims[i,], sens(env=envs[[env]], pars=defpars, par="astar", change=change))
jCPm.sens <- foreach(i=1:n, .combine=rbind) %dopar% with(sims[i,], sens(env=envs[[env]], pars=defpars, par="jCPm", change=change))
jSGm.sens <- foreach(i=1:n, .combine=rbind) %dopar% with(sims[i,], sens(env=envs[[env]], pars=defpars, par="jSGm", change=change))
b.sens <- foreach(i=1:n, .combine=rbind) %dopar% with(sims[i,], sens(env=envs[[env]], pars=defpars, par="b", change=change))

stopCluster(cl)  # Stop cluster

# Combine output with input
sens.all <- ls(pattern="*.sens")
for (i in 1:length(sens.all)) {
  assign(sens.all[i], cbind(sims, get(sens.all[i])))
}

# When growth rate becomes negative, set all responses to NA
for (i in 1:length(sens.all)) {
  if (any(get(sens.all[i])$grchange < 0, na.rm=T)) {
    assign(sens.all[i], within(get(sens.all[i]), {
      shchange <- replace(shchange, grchange < 0, NA)
      cROSchange <- replace(cROSchange, grchange < 0, NA)
      grchange <- replace(grchange, grchange < 0, NA)
    }))
  }
}

# Plot Figure 3
png("img/Fig3.png", width=5, height=5, units="in", res=300)
par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=c(2,0.5,0))
sens.plot(response=c("grchange", "shchange"), pars=c("jCPm", "kCO2", "astar"), cols=c("black", "red", "gold"))
sens.plot(response=c("grchange", "shchange"), pars=c("jNm", "KN"), cols=c("black", "blue"))
sens.plot(response=c("grchange", "shchange"), pars=c("jST0", "jHT0", "jSGm", "jHGm"), cols=c("burlywood4", "burlywood3", "darkolivegreen4", "darkolivegreen3"))
sens.plot(response=c("grchange", "shchange"), pars=c("kNPQ", "kROS", "b"), cols=c("black", "orange", "purple"))
dev.off()