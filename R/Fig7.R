# Load functions
sapply(c("R/def_pars.R",
         "R/init_env.R",
         "R/run_coral.R",
         "R/run_coral_ss.R"), 
       source, .GlobalEnv)

# Increase light for 8000 days then decrease for 8000 days
# Set run time vector
time <- seq(1, 8000, 0.1)  # if single model run, use time input
# Initialize environment
env1 <- init_env(time=time, L=c(25,50,0), N=c(2e-7,2e-7,0), X=c(2e-6,2e-6,0))
L <- list(c(seq(20,40, length.out=40000), seq(40,20, length.out=40000))[-(1:9)])
env1 <- replace(env1, "L", L)

ss <- with(run_coral_ss(env=list(L=20, N=2e-7, X=2e-6), pars=defpars), last(S$S/H$H))

run1 <- run_coral(time=time, env=env1, pars=replace(defpars, "initS", ss))#

png("img/Fig7.png", width=3, height=3, units="in", res=300)
par(mar=c(2,2,1,1), mgp=c(1,0,0), tcl=0.25, cex.axis=0.5, cex.lab=0.75, xaxs="i")
with(run1, {
  plot(env$L, S$S/H$H, type="l", xlab="Light", ylab="S:H biomass", ylim=c(0,0.24))
  for (i in seq(3750,length(time), 1900)) {
    arrows(env$L[i-10], S$S[i-10]/H$H[i-10], 
           env$L[i], S$S[i]/H$H[i], 
           code=2, length=0.075)
  }
})
dev.off()
