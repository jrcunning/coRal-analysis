# Figure 6
# Bleaching in response to increasing light
library(coRal)
library(dplyr)

# Set run time vector
time <- seq(-2, 80, 0.1)  # if single model run, use time input

# Initialize environment
env <- coRal::init_env(time=time, L=c(30,50,0), N=c(1e-7,1e-7,0), X=c(0e-6,0e-6,0))

# Set parameters
pars <- coRal::def_pars()  # Get default parameters

# Run simulation
ss <- with(coRal::run_coral_ss(env=list(L=30, N=1e-7, X=0e-6), pars=pars, dt=0.1), last(S$S/H$H))
run <- coRal::run_coral(time=time, env=env, pars=replace(pars, "initS", ss))

# Plot results
png("img/Fig6.png", width=3, height=3, units="in", res=300)
  with(run, {
    # Set up graphical parameters
    par(mar=c(1.5,1.5,1,1), oma=c(0,0,0,0), mfcol=c(3,2), mgp=c(0.8,0,0), tck=0.025, lwd=1, xaxs="i", 
        cex.main=0.9, cex.axis=0.6, cex.lab=0.7)
    
    # External irradiance
    plot(time, env$L, type="l", col="gold", ylim=c(30,50), xlim=c(0,80), lwd=2, xlab="", ylab=expression(mol~photons~m^{-2}~d^{-1}))
    title("A. Light", adj=0, line=0.25)
    
    # Plot ROS
    plot(NA, xlab="", ylim=c(1, max(2, max(cROS))), xlim=c(0,80), ylab="Relative")
    title(expression(bold("B. ROS production"~italic((c[ROS])))), adj=0, line=0.5)
    lines(time, cROS, col="orange", lwd=2)
    
    # Photosynthesis rate
    plot(NA, xlab="Days", ylim=c(0, max(jCP)), xlim=c(0,80), ylab=expression(mol~C~C-mol~S~d^{-1}))
    title(expression(bold("C. Photosynthesis"~italic((j[CP])))), adj=0, line=0.5)
    lines(time, jCP, col="red", lwd=2)
    
    # Biomass formation - substrate-limitation
    sl <- log(  pmin((rhoN*H/S + rNS)/pars$nNS, pars$jSGm) / pmin(jCP, pars$jSGm) )
    hl <- log(  pmin((jN + pars$nNX*jX + rNH) / pars$nNH, pars$jHGm) / pmin(rhoC*S/H + jX, pars$jHGm) )
    maxabs <- max(c(abs(hl), abs(sl)))
    range <- range(c(hl, sl))
    plot(NA, xlim=c(0,80), xlab="", ylab="Relative", ylim=range)
    title("D. Biomass C-/N-limitation", adj=0, line=0.25)
    lines(time, hl, col="black", lty=1, lwd=1)
    lines(time, sl, col="black", lty=2, lwd=1)
    abline(h=0, lty=3)
    text(par("usr")[2], 0, labels="C-lim.\nN-lim.", cex=0.75, adj=1)
    legend("topleft", legend=c("Host", "Sym"), lty=c(1,2), lwd=c(1,1), col="black", bty="n", cex=0.75)
    
    # Photosynthesis: substrate-limitation
    co2l <- log(  pmin(jL * pars$yCL, pars$jCPm)  / pmin((jCO2 + rCH)*H/S + rCS, pars$jCPm)    )
    maxabs <- max(abs(co2l))
    plot(NA, xlim=c(0,80), xlab="", ylab="Relative", ylim=c(0, maxabs))
    title(expression(bold("E."~CO[2]-limitation)), adj=0, line=0.5)
    lines(time, co2l, col="black", lty=1, lwd=2)
    #abline(h=0, lty=3)
    
    # Symbiont to host biomass ratio
    totSH <- S / H
    totSHf <- totSH[length(totSH)]
    plot(time, totSH, type="l", col="black", ylim=c(0, max(totSH, na.rm=T)), xlim=c(0,80), 
         ylab=expression(C-mol~S~C-mol~H^{-1}), xlab="Days", lwd=2)
    title("F. S:H biomass", adj=0, line=0.25)
  })
dev.off()

