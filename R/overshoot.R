# Load functions
sapply(c("R/def_pars.R",
         "R/init_env.R",
         "R/run_coral.R",
         "R/run_coral_ss.R"), 
       source, .GlobalEnv)

# Load default parameters
defpars <- def_pars()

# Set plotting space for 3x3 panel figure
par(mfcol=c(3,3))

# Run and plot simulations for values of tf (final time/rate of change) and dt (time step)
for (tf in c(500, 1000, 2000)) {
  for (dt in c(0.1, 0.5, 1)) {
    # Set run time vector
    time <- seq(1, tf, dt)  # if single model run, use time input
    
    # Initialize environment (with null light profile - to be replaced)
    env1 <- init_env(time=time, L=c(0,0,0), N=c(1e-7,1e-7,0), X=c(2e-6,2e-6,0))
    
    # Create light profile: increase from 25-50 for 2500 days (=.01 mol/d) then decrease for 2500 days
    L.incr <- seq(20, 50, len=length(time) %/% 2)
    L.decr <- seq(50, 20, len=length(time) - length(L.incr))
    L <- list(c(L.incr, L.decr))
    env1 <- replace(env1, "L", L)
    
    # Calculate rate of change in light
    rate <- diff(range(L))/tf
    
    # Get steady state value of S/H to use as intial value for subsequent simulation
    ss <- with(run_coral_ss(env=list(L=20, N=1e-7, X=2e-6), pars=defpars), last(S$S/H$H))
    
    # Run simulation 
    run1 <- run_coral(time=time, env=env1, pars=replace(defpars, "initS", ss))#
    
    # Plot results
    par(mar=c(2,2,2,1), mgp=c(1,0,0), tcl=0.25, cex.axis=0.5, cex.lab=0.75, xaxs="i", yaxs="i")
    with(run1, {
      plot(env$L, S$S/H$H, type="l", xlab="External light (molph/m2/d)", ylab="S:H biomass", ylim=c(0,0.2))
      text(par("usr")[2], par("usr")[4], adj=c(1,1),
           labels=paste("Rate of change in light:", rate, "/d \n dt:", dt, "d"), cex.main=0.9)
      for (i in seq(100, length(time), last(time)/2)) {
        arrows(env$L[i-1], S$S[i-1]/H$H[i-1], 
               env$L[i], S$S[i]/H$H[i], 
               code=2, length=0.075)
      }
    })
  }
}

