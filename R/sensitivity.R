library(parallel)
library(doParallel)
library(foreach)
source("R/def_pars.R")
source("R/run_coral.R")
source("R/init_env.R")

# Sensitivity analysis

# Create function to measure sensitivity of responses to specified change in single parameter
sens <- function(time, env, pars, par, change) {
  # Run 1 
  run1 <- run_coral(time=time, env=env, pars=pars)
  run1ss <- lapply(run1[c("H", "S")], function(x) x[length(time), ])
  run1gr <- run1ss$H$dH.Hdt
  run1sh <- run1ss$S$S/run1ss$H$H
  # Run 2 - with changed parameter
  run2 <- run_coral(time=time, env=env, pars=replace(pars, par, with(pars, get(par))*change))
  run2ss <- lapply(run2[c("H", "S")], function(x) x[length(time), ])
  run2gr <- run2ss$H$dH.Hdt
  run2sh <- run2ss$S$S/run2ss$H$H
  # Calculate relative change in steady state response
  grchange <- run2gr / run1gr
  shchange <- run2sh / run1sh
  cROSchange <- run2ss$S$cROS / run1ss$S$cROS
  # Return relative changes
  return(list(grchange=grchange, shchange=shchange, cROSchange=cROSchange))
}

# Define function to measure sensitivity across a range of changes in parameter value
vsens <- function(par, levs, env) {
  gr <- rep(NA, length(levs))
  sh <- rep(NA, length(levs))
  cROS <- rep(NA, length(levs))
  for (i in 1:length(levs)) {
    sa <- sens(time=time, env=env, pars=defpars, par=par, change=levs[i])
    gr[i] <- sa[[1]]
    sh[i] <- sa[[2]]
    cROS[i] <- sa[[3]]
  }
  return(data.frame(lev=levs, gr=gr, sh=sh, cROS=cROS))
}

# group parameters by process, then evaluate sensitivity under diff environments and plot together.

# Set up time vector and parameters
time <- seq(1,50,0.5)  # Set time
defpars <- def_pars()  # Get default parameters
# Define four environments for conducting sensitivity analysis
envs <- list(
  LLLN=init_env(time=time, L=c(2,2,0), N=c(1e-7,1e-7,0), X=c(0,0,0)),
  LLHN=init_env(time=time, L=c(2,2,0), N=c(4e-6,4e-6,0), X=c(0,0,0)),
  HLLN=init_env(time=time, L=c(25,25,0), N=c(1e-7,1e-7,0), X=c(0,0,0)),
  HLHN=init_env(time=time, L=c(25,25,0), N=c(4e-6,4e-6,0), X=c(0,0,0))
)

# Set up cluster for parallel processing
cl <- makeCluster(detectCores())  # Initiate cluster
registerDoParallel(cl)

# Run sensitivity analysis in each environment for each parameter
levs <- seq(0.5, 1.5, 0.1)  # Set levels of % change in parameter value to evaluate
jCPm.sens <- foreach(i=1:4) %dopar% vsens(par="jCPm", levs=levs, env=envs[[i]])
astar.sens <- foreach(i=1:4) %dopar% vsens(par="astar", levs=levs, env=envs[[i]])
jNm.sens <- foreach(i=1:4) %dopar% vsens(par="jNm", levs=levs, env=envs[[i]])
KN.sens <- foreach(i=1:4) %dopar% vsens(par="KN", levs=levs, env=envs[[i]])
jST0.sens <- foreach(i=1:4) %dopar% vsens(par="jST0", levs=levs, env=envs[[i]])
jHT0.sens <- foreach(i=1:4) %dopar% vsens(par="jHT0", levs=levs, env=envs[[i]])
jSGm.sens <- foreach(i=1:4) %dopar% vsens(par="jSGm", levs=levs, env=envs[[i]])
jHGm.sens <- foreach(i=1:4) %dopar% vsens(par="jHGm", levs=levs, env=envs[[i]])
kNPQ.sens <- foreach(i=1:4) %dopar% vsens(par="kNPQ", levs=levs, env=envs[[i]])
kROS.sens <- foreach(i=1:4) %dopar% vsens(par="kROS", levs=levs, env=envs[[i]])
kCO2.sens <- foreach(i=1:4) %dopar% vsens(par="kCO2", levs=levs, env=envs[[i]])
b.sens <- foreach(i=1:4) %dopar% vsens(par="b", levs=levs, env=envs[[i]])
#jXm
#KX

stopCluster(cl)  # Stop cluster

names(jCPm.sens) <- names(envs)
names(jCO2a.sens) <- names(envs)
names(astar.sens) <- names(envs)
names(jNm.sens) <- names(envs)
names(KN.sens) <- names(envs)
names(jST0.sens) <- names(envs)
names(jHT0.sens) <- names(envs)
names(jSGm.sens) <- names(envs)
names(jHGm.sens) <- names(envs)
names(kNPQ.sens) <- names(envs)
names(kROS.sens) <- names(envs)


# PLOT
#png("img/Fig3.png", width=5, heigh=5, units="in", res=300)
par(mfrow=c(2,2), mar=c(2,2,2,1), mgp=c(1,0,0), tcl=0.2)
# plot photosynthesis parameters: jCPm, jCO2a, astar
plot(NA, xlim=range(levs), ylim=range(levs),
     xlab="Proportional change in parameter", ylab="Proportional change in growth")
title("A. Photosynthesis parameters", adj=0)
legend("topleft", legend=names(envs), lwd=2, lty=c(1,2,3,4), cex=0.7, bty="n")
legend("bottomright", legend=c("jCPm", "jCO2a", "astar"), cex=0.7, lwd=2, lty=1, col=c("black", "red", "gold"), bty="n")
with(jCPm.sens, {
  lines(LLLN$gr ~ LLLN$lev, lwd=2, lty=1)
  lines(LLHN$gr ~ LLHN$lev, lwd=2, lty=2)
  lines(HLLN$gr ~ HLLN$lev, lwd=2, lty=3)
  lines(HLHN$gr ~ HLHN$lev, lwd=2, lty=4)
})
with(jCO2a.sens, {
  lines(LLLN$gr ~ LLLN$lev, lwd=2, lty=1, col="red")
  lines(LLHN$gr ~ LLHN$lev, lwd=2, lty=2, col="red")
  lines(HLLN$gr ~ HLLN$lev, lwd=2, lty=3, col="red")
  lines(HLHN$gr ~ HLHN$lev, lwd=2, lty=4, col="red")
})
with(astar.sens, {
  lines(LLLN$gr ~ LLLN$lev, lwd=2, lty=1, col="gold")
  lines(LLHN$gr ~ LLHN$lev, lwd=2, lty=2, col="gold")
  lines(HLLN$gr ~ HLLN$lev, lwd=2, lty=3, col="gold")
  lines(HLHN$gr ~ HLHN$lev, lwd=2, lty=4, col="gold")
})

# plot nitrogen parameters: jNm, KN
plot(NA, xlim=range(levs), ylim=range(levs),
     xlab="Proportional change in parameter", ylab="Proportional change in growth")
title("B. DIN uptake parameters", adj=0)
legend("topleft", legend=names(envs), lwd=2, lty=c(1,2,3,4), cex=0.7, bty="n")
legend("bottomright", legend=c("jNm", "KN"), cex=0.7, lwd=2, lty=1, col=c("black", "blue"), bty="n")
with(jNm.sens, {
  lines(LLLN$gr ~ LLLN$lev, lwd=2, lty=1)
  lines(LLHN$gr ~ LLHN$lev, lwd=2, lty=2)
  lines(HLLN$gr ~ HLLN$lev, lwd=2, lty=3)
  lines(HLHN$gr ~ HLHN$lev, lwd=2, lty=4)
})
with(KN.sens, {
  lines(LLLN$gr ~ LLLN$lev, lwd=2, lty=1, col="blue")
  lines(LLHN$gr ~ LLHN$lev, lwd=2, lty=2, col="blue")
  lines(HLLN$gr ~ HLLN$lev, lwd=2, lty=3, col="blue")
  lines(HLHN$gr ~ HLHN$lev, lwd=2, lty=4, col="blue")
})

# plot biomass parameters: jST0, jHT0, jSGm, jHGm
plot(NA, xlim=range(levs), ylim=range(levs),
     xlab="Proportional change in parameter", ylab="Proportional change in growth")
title("C. Biomass parameters", adj=0)
legend("topleft", legend=names(envs), lwd=2, lty=c(1,2,3,4), cex=0.7, bty="n")
legend("bottomright", legend=c("jHGm","jHT0","jSGm","jST0"), cex=0.7, lwd=2, lty=1, col=c("burlywood4", "burlywood3", "darkolivegreen4", "darkolivegreen3"), bty="n")
with(jHGm.sens, {
  lines(LLLN$gr ~ LLLN$lev, lwd=2, lty=1, col="burlywood4")
  lines(LLHN$gr ~ LLHN$lev, lwd=2, lty=2, col="burlywood4")
  lines(HLLN$gr ~ HLLN$lev, lwd=2, lty=3, col="burlywood4")
  lines(HLHN$gr ~ HLHN$lev, lwd=2, lty=4, col="burlywood4")
})
with(jHT0.sens, {
  lines(LLLN$gr ~ LLLN$lev, lwd=2, lty=1, col="burlywood3")
  lines(LLHN$gr ~ LLHN$lev, lwd=2, lty=2, col="burlywood3")
  lines(HLLN$gr ~ HLLN$lev, lwd=2, lty=3, col="burlywood3")
  lines(HLHN$gr ~ HLHN$lev, lwd=2, lty=4, col="burlywood3")
})
with(jSGm.sens, {
  lines(LLLN$gr ~ LLLN$lev, lwd=2, lty=1, col="darkolivegreen4")
  lines(LLHN$gr ~ LLHN$lev, lwd=2, lty=2, col="darkolivegreen4")
  lines(HLLN$gr ~ HLLN$lev, lwd=2, lty=3, col="darkolivegreen4")
  lines(HLHN$gr ~ HLHN$lev, lwd=2, lty=4, col="darkolivegreen4")
})
with(jST0.sens, {
  lines(LLLN$gr ~ LLLN$lev, lwd=2, lty=1, col="darkolivegreen3")
  lines(LLHN$gr ~ LLHN$lev, lwd=2, lty=2, col="darkolivegreen3")
  lines(HLLN$gr ~ HLLN$lev, lwd=2, lty=3, col="darkolivegreen3")
  lines(HLHN$gr ~ HLHN$lev, lwd=2, lty=4, col="darkolivegreen3")
})

# Plot stress parameters: kNPQ, kROS
plot(NA, xlim=range(levs), ylim=range(levs),
     xlab="Proportional change in parameter", ylab="Proportional change in growth")
title("D. Stress parameters", adj=0)
legend("topleft", legend=names(envs), lwd=2, lty=c(1,2,3,4), cex=0.7, bty="n")
legend("bottomright", legend=c("kNPQ", "kROS"), cex=0.7, lwd=2, lty=1, col=c("black", "blue"), bty="n")
with(kNPQ.sens, {
  lines(LLLN$gr ~ LLLN$lev, lwd=2, lty=1)
  lines(LLHN$gr ~ LLHN$lev, lwd=2, lty=2)
  lines(HLLN$gr ~ HLLN$lev, lwd=2, lty=3)
  lines(HLHN$gr ~ HLHN$lev, lwd=2, lty=4)
})
with(kROS.sens, {
  lines(LLLN$gr ~ LLLN$lev, lwd=2, lty=1, col="blue")
  lines(LLHN$gr ~ LLHN$lev, lwd=2, lty=2, col="blue")
  lines(HLLN$gr ~ HLLN$lev, lwd=2, lty=3, col="blue")
  lines(HLHN$gr ~ HLHN$lev, lwd=2, lty=4, col="blue")
})
#dev.off()


# PLOT

par(mfrow=c(2,2), mar=c(2,2,2,1), mgp=c(1,0,0), tcl=0.2)
# plot photosynthesis parameters: jCPm, jCO2a, astar
plot(NA, xlim=range(levs), ylim=range(levs),
     xlab="Proportional change in parameter", ylab="Proportional change in S:H")
title("A. Photosynthesis parameters", adj=0)
legend("topleft", legend=names(envs), lwd=2, lty=c(1,2,3,4), cex=0.7, bty="n")
legend("bottomright", legend=c("jCPm", "jCO2a", "astar"), cex=0.7, lwd=2, lty=1, col=c("black", "red", "gold"), bty="n")
with(jCPm.sens, {
  lines(LLLN$sh ~ LLLN$lev, lwd=2, lty=1)
  lines(LLHN$sh ~ LLHN$lev, lwd=2, lty=2)
  lines(HLLN$sh ~ HLLN$lev, lwd=2, lty=3)
  lines(HLHN$sh ~ HLHN$lev, lwd=2, lty=4)
})
with(jCO2a.sens, {
  lines(LLLN$sh ~ LLLN$lev, lwd=2, lty=1, col="red")
  lines(LLHN$sh ~ LLHN$lev, lwd=2, lty=2, col="red")
  lines(HLLN$sh ~ HLLN$lev, lwd=2, lty=3, col="red")
  lines(HLHN$sh ~ HLHN$lev, lwd=2, lty=4, col="red")
})
with(astar.sens, {
  lines(LLLN$sh ~ LLLN$lev, lwd=2, lty=1, col="gold")
  lines(LLHN$sh ~ LLHN$lev, lwd=2, lty=2, col="gold")
  lines(HLLN$sh ~ HLLN$lev, lwd=2, lty=3, col="gold")
  lines(HLHN$sh ~ HLHN$lev, lwd=2, lty=4, col="gold")
})

# plot nitrogen parameters: jNm, KN
plot(NA, xlim=range(levs), ylim=range(levs),
     xlab="Proportional change in parameter", ylab="Proportional change in S:H")
title("B. DIN uptake parameters", adj=0)
legend("topleft", legend=names(envs), lwd=2, lty=c(1,2,3,4), cex=0.7, bty="n")
legend("bottomright", legend=c("jNm", "KN"), cex=0.7, lwd=2, lty=1, col=c("black", "blue"), bty="n")
with(jNm.sens, {
  lines(LLLN$sh ~ LLLN$lev, lwd=2, lty=1)
  lines(LLHN$sh ~ LLHN$lev, lwd=2, lty=2)
  lines(HLLN$sh ~ HLLN$lev, lwd=2, lty=3)
  lines(HLHN$sh ~ HLHN$lev, lwd=2, lty=4)
})
with(KN.sens, {
  lines(LLLN$sh ~ LLLN$lev, lwd=2, lty=1, col="blue")
  lines(LLHN$sh ~ LLHN$lev, lwd=2, lty=2, col="blue")
  lines(HLLN$sh ~ HLLN$lev, lwd=2, lty=3, col="blue")
  lines(HLHN$sh ~ HLHN$lev, lwd=2, lty=4, col="blue")
})

# plot biomass parameters: jST0, jHT0, jSGm, jHGm
plot(NA, xlim=range(levs), ylim=range(levs),
     xlab="Proportional change in parameter", ylab="Proportional change in S:H")
title("C. Biomass parameters", adj=0)
legend("topleft", legend=names(envs), lwd=2, lty=c(1,2,3,4), cex=0.7, bty="n")
legend("bottomright", legend=c("jHGm","jHT0","jSGm","jST0"), cex=0.7, lwd=2, lty=1, col=c("burlywood4", "burlywood3", "darkolivegreen4", "darkolivegreen3"), bty="n")
with(jHGm.sens, {
  lines(LLLN$sh ~ LLLN$lev, lwd=2, lty=1, col="burlywood4")
  lines(LLHN$sh ~ LLHN$lev, lwd=2, lty=2, col="burlywood4")
  lines(HLLN$sh ~ HLLN$lev, lwd=2, lty=3, col="burlywood4")
  lines(HLHN$sh ~ HLHN$lev, lwd=2, lty=4, col="burlywood4")
})
with(jHT0.sens, {
  lines(LLLN$sh ~ LLLN$lev, lwd=2, lty=1, col="burlywood3")
  lines(LLHN$sh ~ LLHN$lev, lwd=2, lty=2, col="burlywood3")
  lines(HLLN$sh ~ HLLN$lev, lwd=2, lty=3, col="burlywood3")
  lines(HLHN$sh ~ HLHN$lev, lwd=2, lty=4, col="burlywood3")
})
with(jSGm.sens, {
  lines(LLLN$sh ~ LLLN$lev, lwd=2, lty=1, col="darkolivegreen4")
  lines(LLHN$sh ~ LLHN$lev, lwd=2, lty=2, col="darkolivegreen4")
  lines(HLLN$sh ~ HLLN$lev, lwd=2, lty=3, col="darkolivegreen4")
  lines(HLHN$sh ~ HLHN$lev, lwd=2, lty=4, col="darkolivegreen4")
})
with(jST0.sens, {
  lines(LLLN$sh ~ LLLN$lev, lwd=2, lty=1, col="darkolivegreen3")
  lines(LLHN$sh ~ LLHN$lev, lwd=2, lty=2, col="darkolivegreen3")
  lines(HLLN$sh ~ HLLN$lev, lwd=2, lty=3, col="darkolivegreen3")
  lines(HLHN$sh ~ HLHN$lev, lwd=2, lty=4, col="darkolivegreen3")
})

# Plot stress parameters: kNPQ, kROS
plot(NA, xlim=range(levs), ylim=range(levs),
     xlab="Proportional change in parameter", ylab="Proportional change in S:H")
title("D. Stress parameters", adj=0)
legend("topleft", legend=names(envs), lwd=2, lty=c(1,2,3,4), cex=0.7, bty="n")
legend("bottomright", legend=c("kNPQ", "kROS"), cex=0.7, lwd=2, lty=1, col=c("black", "blue"), bty="n")
with(kNPQ.sens, {
  lines(LLLN$sh ~ LLLN$lev, lwd=2, lty=1)
  lines(LLHN$sh ~ LLHN$lev, lwd=2, lty=2)
  lines(HLLN$sh ~ HLLN$lev, lwd=2, lty=3)
  lines(HLHN$sh ~ HLHN$lev, lwd=2, lty=4)
})
with(kROS.sens, {
  lines(LLLN$sh ~ LLLN$lev, lwd=2, lty=1, col="blue")
  lines(LLHN$sh ~ LLHN$lev, lwd=2, lty=2, col="blue")
  lines(HLLN$sh ~ HLLN$lev, lwd=2, lty=3, col="blue")
  lines(HLHN$sh ~ HLHN$lev, lwd=2, lty=4, col="blue")
})




# PLOT CROS
# PLOT
#png("img/Fig3.png", width=5, heigh=5, units="in", res=300)
par(mfrow=c(2,2), mar=c(2,2,2,1), mgp=c(1,0,0), tcl=0.2)
# plot photosynthesis parameters: jCPm, jCO2a, astar
plot(NA, xlim=range(levs), ylim=range(levs),
     xlab="Proportional change in parameter", ylab="Proportional change in growth")
title("A. Photosynthesis parameters", adj=0)
legend("topleft", legend=names(envs), lwd=2, lty=c(1,2,3,4), cex=0.7, bty="n")
legend("bottomright", legend=c("jCPm", "jCO2a", "astar"), cex=0.7, lwd=2, lty=1, col=c("black", "red", "gold"), bty="n")
with(jCPm.sens, {
  lines(LLLN$cROS ~ LLLN$lev, lwd=2, lty=1)
  lines(LLHN$cROS ~ LLHN$lev, lwd=2, lty=2)
  lines(HLLN$cROS ~ HLLN$lev, lwd=2, lty=3)
  lines(HLHN$cROS ~ HLHN$lev, lwd=2, lty=4)
})
with(jCO2a.sens, {
  lines(LLLN$cROS ~ LLLN$lev, lwd=2, lty=1, col="red")
  lines(LLHN$cROS ~ LLHN$lev, lwd=2, lty=2, col="red")
  lines(HLLN$cROS ~ HLLN$lev, lwd=2, lty=3, col="red")
  lines(HLHN$cROS ~ HLHN$lev, lwd=2, lty=4, col="red")
})
with(astar.sens, {
  lines(LLLN$cROS ~ LLLN$lev, lwd=2, lty=1, col="gold")
  lines(LLHN$cROS ~ LLHN$lev, lwd=2, lty=2, col="gold")
  lines(HLLN$cROS ~ HLLN$lev, lwd=2, lty=3, col="gold")
  lines(HLHN$cROS ~ HLHN$lev, lwd=2, lty=4, col="gold")
})

# plot nitrogen parameters: jNm, KN
plot(NA, xlim=range(levs), ylim=range(levs),
     xlab="Proportional change in parameter", ylab="Proportional change in growth")
title("B. DIN uptake parameters", adj=0)
legend("topleft", legend=names(envs), lwd=2, lty=c(1,2,3,4), cex=0.7, bty="n")
legend("bottomright", legend=c("jNm", "KN"), cex=0.7, lwd=2, lty=1, col=c("black", "blue"), bty="n")
with(jNm.sens, {
  lines(LLLN$cROS ~ LLLN$lev, lwd=2, lty=1)
  lines(LLHN$cROS ~ LLHN$lev, lwd=2, lty=2)
  lines(HLLN$cROS ~ HLLN$lev, lwd=2, lty=3)
  lines(HLHN$cROS ~ HLHN$lev, lwd=2, lty=4)
})
with(KN.sens, {
  lines(LLLN$cROS ~ LLLN$lev, lwd=2, lty=1, col="blue")
  lines(LLHN$cROS ~ LLHN$lev, lwd=2, lty=2, col="blue")
  lines(HLLN$cROS ~ HLLN$lev, lwd=2, lty=3, col="blue")
  lines(HLHN$cROS ~ HLHN$lev, lwd=2, lty=4, col="blue")
})

# plot biomass parameters: jST0, jHT0, jSGm, jHGm
plot(NA, xlim=range(levs), ylim=range(levs),
     xlab="Proportional change in parameter", ylab="Proportional change in growth")
title("C. Biomass parameters", adj=0)
legend("topleft", legend=names(envs), lwd=2, lty=c(1,2,3,4), cex=0.7, bty="n")
legend("bottomright", legend=c("jHGm","jHT0","jSGm","jST0"), cex=0.7, lwd=2, lty=1, col=c("burlywood4", "burlywood3", "darkolivegreen4", "darkolivegreen3"), bty="n")
with(jHGm.sens, {
  lines(LLLN$cROS ~ LLLN$lev, lwd=2, lty=1, col="burlywood4")
  lines(LLHN$cROS ~ LLHN$lev, lwd=2, lty=2, col="burlywood4")
  lines(HLLN$cROS ~ HLLN$lev, lwd=2, lty=3, col="burlywood4")
  lines(HLHN$cROS ~ HLHN$lev, lwd=2, lty=4, col="burlywood4")
})
with(jHT0.sens, {
  lines(LLLN$cROS ~ LLLN$lev, lwd=2, lty=1, col="burlywood3")
  lines(LLHN$cROS ~ LLHN$lev, lwd=2, lty=2, col="burlywood3")
  lines(HLLN$cROS ~ HLLN$lev, lwd=2, lty=3, col="burlywood3")
  lines(HLHN$cROS ~ HLHN$lev, lwd=2, lty=4, col="burlywood3")
})
with(jSGm.sens, {
  lines(LLLN$cROS ~ LLLN$lev, lwd=2, lty=1, col="darkolivegreen4")
  lines(LLHN$cROS ~ LLHN$lev, lwd=2, lty=2, col="darkolivegreen4")
  lines(HLLN$cROS ~ HLLN$lev, lwd=2, lty=3, col="darkolivegreen4")
  lines(HLHN$cROS ~ HLHN$lev, lwd=2, lty=4, col="darkolivegreen4")
})
with(jST0.sens, {
  lines(LLLN$cROS ~ LLLN$lev, lwd=2, lty=1, col="darkolivegreen3")
  lines(LLHN$cROS ~ LLHN$lev, lwd=2, lty=2, col="darkolivegreen3")
  lines(HLLN$cROS ~ HLLN$lev, lwd=2, lty=3, col="darkolivegreen3")
  lines(HLHN$cROS ~ HLHN$lev, lwd=2, lty=4, col="darkolivegreen3")
})

# Plot stress parameters: kNPQ, kROS
plot(NA, xlim=range(levs), ylim=range(levs),
     xlab="Proportional change in parameter", ylab="Proportional change in growth")
title("D. Stress parameters", adj=0)
legend("topleft", legend=names(envs), lwd=2, lty=c(1,2,3,4), cex=0.7, bty="n")
legend("bottomright", legend=c("kNPQ", "kROS"), cex=0.7, lwd=2, lty=1, col=c("black", "blue"), bty="n")
with(kNPQ.sens, {
  lines(LLLN$cROS ~ LLLN$lev, lwd=2, lty=1)
  lines(LLHN$cROS ~ LLHN$lev, lwd=2, lty=2)
  lines(HLLN$cROS ~ HLLN$lev, lwd=2, lty=3)
  lines(HLHN$cROS ~ HLHN$lev, lwd=2, lty=4)
})
with(kROS.sens, {
  lines(LLLN$cROS ~ LLLN$lev, lwd=2, lty=1, col="blue")
  lines(LLHN$cROS ~ LLHN$lev, lwd=2, lty=2, col="blue")
  lines(HLLN$cROS ~ HLLN$lev, lwd=2, lty=3, col="blue")
  lines(HLHN$cROS ~ HLHN$lev, lwd=2, lty=4, col="blue")
})

