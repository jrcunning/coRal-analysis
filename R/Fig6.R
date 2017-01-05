# Figure 5
# Bleaching in response to increasing light

# Load functions
sapply(c("R/def_pars.R",
         "R/init_env.R",
         "R/run_coral.R",
         "R/run_coral_ss.R"), 
       source, .GlobalEnv)

# Set run time vector
time <- seq(1, 80, 0.1)  # if single model run, use time input

# Initialize environment
env <- init_env(time=time, L=c(30,50,0), N=c(1e-7,1e-7,0), X=c(0e-6,0e-6,0))
#L <- list(c(rep(20,100), rep(50,300))[-(1:9)])
#env <- replace(env, "L", L)

# Set parameters
defpars <- def_pars()  # Get default parameters

# Run simulation
ss <- with(run_coral_ss(env=list(L=25, N=1e-7, X=0e-6), pars=defpars, dt=0.1), last(S$S/H$H))
run <- run_coral(time=time, env=env, pars=replace(defpars, "initS", ss))

# Plot results
png("img/Fig6.png", width=3, height=3, units="in", res=300)
  with(run, {
    # Set up graphical parameters
    par(mar=c(1.5,1.5,1,1), oma=c(0,0,0,0), mfcol=c(3,2), mgp=c(0.8,0,0), tck=0.025, lwd=1, xaxs="i", 
        cex.main=0.9, cex.axis=0.6, cex.lab=0.7)
    
    # External irradiance
    plot(time, env$L, type="l", col="gold", ylim=c(30,50), lwd=2, xlab="", ylab="molph/m2/d")
    title("A. External irradiance", adj=0, line=0.25)
    
    # Plot ROS
    plot(NA, xlim=range(time), xlab="", ylim=c(1, max(2, max(S$cROS))), ylab="Relative")
    title("B. ROS production", adj=0.05, line=0)
    lines(time, S$cROS, col="orange", lwd=2)
    
    # Photosynthesis rate
    plot(NA, xlim=range(time), xlab="Days", ylim=c(0, max(S$jCP)), ylab="molC/CmolS/d")
    title("C. Photosynthesis rate", adj=0.05, line=0)
    lines(time, S$jCP, col="red", lwd=2)
    
    # Biomass formation - substrate-limitation
    sl <- log(  pmin((H$rhoN*H$H/S$S + S$rNS)/pars$nNS, pars$jSGm) / pmin(S$jCP, pars$jSGm) )
    hl <- log(  pmin((H$jN + pars$nNX*H$jX + H$rNH) / pars$nNH, pars$jHGm) / pmin(S$rhoC*S$S/H$H + H$jX, pars$jHGm) )
    maxabs <- max(c(abs(hl), abs(sl)))
    range <- range(c(hl, sl))
    plot(NA, xlim=range(time), xlab="", ylab="Relative", ylim=range)
    title("D. Biomass C-/N-limitation", adj=0, line=0.25)
    lines(time, hl, col="black", lty=1, lwd=1)
    lines(time, sl, col="black", lty=2, lwd=1)
    abline(h=0, lty=3)
    text(par("usr")[2], 0, labels="C-lim.\nN-lim.", cex=0.75, adj=1)
    legend("topleft", legend=c("Host", "Sym"), lty=c(1,2), lwd=c(1,1), col="black", bty="n", cex=0.75)
    
    # Photosynthesis: substrate-limitation
    co2l <- log(  pmin(S$jL * pars$yCL, pars$jCPm)  / pmin((H$jCO2 + H$rCH)*H$H/S$S + S$rCS, pars$jCPm)    )
    maxabs <- max(abs(co2l))
    plot(NA, xlim=range(time), xlab="", ylab="Relative", ylim=c(0, maxabs))
    title("E. CO2-limitation", adj=0, line=0.25)
    lines(time, co2l, col="black", lty=1, lwd=2)
    #abline(h=0, lty=3)
    
    # Symbiont to host biomass ratio
    totSH <- S$S / H$H
    totSHf <- totSH[length(totSH)]
    plot(time, totSH, type="l", col="black", ylim=c(0, max(totSH, na.rm=T)), ylab="CmolS/CmolH", xlab="Days", lwd=2)
    title("F. S:H biomass", adj=0, line=0.25)
  })
dev.off()

