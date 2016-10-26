# Load functions
sapply(c("R/def_pars.R",
         "R/init_env.R",
         "R/run_coral.R",
         "R/run_coral_ss.R"), 
       source, .GlobalEnv)

# Load default parameters
defpars <- def_pars()

# Set run time vector
time <- seq(1, 5000, 0.1)  # if single model run, use time input

# Initialize environment (with null light profile - to be replaced)
env1 <- init_env(time=time, L=c(0,0,0), N=c(2e-7,2e-7,0), X=c(2e-6,2e-6,0))

# Create light profile: increase from 20-40 for 2500 days (=.008 mol/d) then decrease for 2500 days
L.incr <- seq(20, 40, len=length(time) %/% 2)
L.decr <- seq(40, 20, len=length(time) - length(L.incr))
L <- list(c(L.incr, L.decr))
env1 <- replace(env1, "L", L)

ss <- with(run_coral_ss(env=list(L=20, N=2e-7, X=2e-6), pars=defpars), last(S$S/H$H))

run1 <- run_coral(time=time, env=env1, pars=replace(defpars, "initS", ss))#

png("img/Fig7.png", width=4, height=4, units="in", res=300)
par(mar=c(2,2,1,1), mgp=c(1,0,0), tcl=0.25, cex.axis=0.5, cex.lab=0.75, xaxs="i")
with(run1, {
  plot(env$L, S$S/H$H, type="l", xlab="Light", ylab="S:H biomass", ylim=c(0,0.24))
  for (i in seq(1300, length(time), last(time)/3.5)) {
    arrows(env$L[i-10], S$S[i-10]/H$H[i-10], 
           env$L[i], S$S[i]/H$H[i], 
           code=2, length=0.075)
  }
})
dev.off()
