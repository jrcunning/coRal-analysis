# This script contains functions that can be used to plot the results of a model run.
# Fluxes can be plotted individually, or altogether with the "plot_run" function.
# The argument to these functions should be an object returned from "run_coral.R".

# Plotting functions for individual fluxes
# ========================================

# External irradiance
plot_L <- function(run) with(run, {
  plot(time, env$L, type="l", col="yellow", ylim=c(0,60), lwd=3, xlab="", ylab="mol/m2/d", main="Irradiance")
})

# External DIN concentration
plot_DIN <- function(run) with(run, {
  plot(time, env$N * 10^6, type="l", col="blue", ylim=c(0,4), lwd=3, xlab="", ylab="µmol/L", main="DIN")
})

# Specific growth rates of host and symbiont
plot_gr <- function(run) with(run, {
  Hgrf <- H$dH.Hdt[length(H$dH.Hdt)]
  plot(time, H$dH.Hdt, type="l", ylim=c(min(0, min(c(H$dH.Hdt[-(1:10)], S$dS.Sdt[-(1:10)]))), max(c(H$dH.Hdt[-(1:10)], S$dS.Sdt[-(1:10)]))),
       xlab="", ylab="d-1", lwd=3, cex=1, cex.lab=1, 
       main="Specific growth rate")
  if(any(c(H$dH.Hdt[-(1:10)], S$dS.Sdt[-(1:10)])<0)) abline(h=0, lwd=0.5, lty=3)
  text(time[0.95*length(time)], Hgrf, labels=as.character(round(Hgrf, 3)), pos=3, xpd=T, cex=0.75)
  lines(time, S$dS.Sdt, col="black", lwd=1)
  legend("topright", legend=c("Host", "Sym"), lwd=c(3,1), col="black", bty="n", y.intersp=1)
})

# Symbiont to host ratio
plot_sh <- function(run) with(run, {
  totSH <- S$S / H$H
  totSHf <- totSH[length(totSH)]
  plot(time, totSH, type="l", col="black", ylim=c(0, max(totSH, na.rm=T)), ylab="C-molS/C-molH", xlab="", lwd=3,
       main="Symbiont:host biomass")
  text(time[0.95*length(time)], totSHf, labels=as.character(round(totSHf, 3)), pos=3, xpd=T, cex=0.75)
})

# Light quenching (carbon fixation, NPQ, and excess (=ROS producing))
plot_Lq <- function(run) with(run, {
  plot(NA, xlim=range(time), xlab="", ylim=c(0, max(c(S$jL, S$jeL))), ylab="mol photons/C-molS/d",
       main="Light quenching")
  lines(time, S$jL, col="yellow", lty=3, lwd=3) # total amount absorbed
  lines(time, S$jL - S$jeL, col="yellow", lwd=1) # amt. used in photosynthesis
  lines(time, (S$jL - S$jeL) + pmin(S$jeL, pars$jNPQ), col="yellow", lty=1, lwd=3) # amt. quenched by NPQ
  legend("topright", legend=c("Excess", "NPQ", "Photo."), lty=c(3,1,1), lwd=c(3,3,1), col="yellow", bty="n")
})

# Photosynthesis rate (carbon fixation)
plot_photo <- function(run) with(run, {
  plot(NA, xlim=range(time), xlab="", ylim=c(0, max(S$jCP)), ylab="C-mol/C-molS/d")
  title("Photosynthesis rate", adj=0.05, line=0)
  lines(time, S$jCP, col="red", lwd=3)
})

# Photosynthesis SU dynamics (proportions of light and CO2 rejected from SU)
plot_PSU <- function(run) with(run, {
  plot(NA, xlim=range(time), xlab="", ylim=c(0,1), ylab="Proportion of uptake")
  title("Photosynthesis: substrate excess", adj=0.05, line=0)
  lines(time, S$jeL/S$jL, col="yellow", lwd=3)
  lines(time, S$jCO2w/((H$jCO2 + H$rCH)*H$H/S$S + S$rCS), col="red", lwd=3)
  legend("topright", legend=c("Light", "DIC"), lty=c(1,1), lwd=3, bty="n", col=c("yellow", "red"))
})

# Relative ROS production due to excitation energy in excess of carbon fixation and NPQ
plot_ROS <- function(run) with(run, {
  plot(NA, xlim=range(time), xlab="", ylim=c(1, max(S$cROS)), ylab="Relative to baseline")
  title("ROS production", adj=0.05, line=0)
  lines(time, S$cROS, col="orange", lwd=3)
})

# Symbiont biomass SU dynamics (proportions of C and N rejected from SU)
plot_symSU <- function(run) with(run, {
  plot(NA, xlim=range(time), xlab="", ylim=c(0, 1), ylab="Proportion of uptake")
  title("Symbiont: substrate excess", adj=0.05, line=0)
  lines(time, S$rhoC/S$jCP, col="red", lwd=3)
  lines(time, S$jNw/(H$rhoN*H$H/S$S + S$rNS), col="blue", lwd=3)
  legend("topright", legend=c("C", "N"), lwd=3, bty="n", col=c("red", "blue"))
})

# Host biomass SU dynamics (proportions of C and N rejected from SU)
plot_corSU <- function(run) with(run, {
  plot(NA, xlim=range(time), xlab="", ylim=c(0, 1), ylab="Proportion of uptake")
  title("Coral: substrate excess", adj=0.05, line=0)
  lines(time, H$jCw/(S$rhoC*S$S/H$H + H$jX), col="red", lwd=3)
  lines(time, H$rhoN/(H$jN + H$jX*pars$nNX + H$rNH), col="blue", lwd=3)
  legend("topright", legend=c("C", "N"), lwd=3, bty="n", col=c("red", "blue"))
})



# Plotting function for all fluxes (as multi-panel figure)
# ========================================================

plot_run <- function(run) {
  with(run, {
    # Set up graphical parameters
    par(mar=c(2,2,1,1), oma=c(0,0,0,0), mfrow=c(5,2), mgp=c(1.2,0.2,0), tck=0.025, lwd=1, xaxs="i", 
        cex.main=1.5, cex.axis=1, cex.lab=1)
    
    # External irradiance
    plot(time, env$L, type="l", col="yellow", ylim=c(0,60), lwd=3, xlab="", ylab="mol/m2/d")
    title("Light", adj=0.05, line=0)
    
    # External DIN concentration
    plot(time, env$N * 10^6, type="l", col="blue", ylim=c(0,4), ylab="µmol/L", xlab="", lwd=3)
    title("DIN", adj=0.05, line=0)
    
    # Specific growth rates of host and symbiont
    Hgrf <- H$dH.Hdt[length(H$dH.Hdt)]
    plot(time, H$dH.Hdt, type="l", ylim=c(min(0, min(c(H$dH.Hdt[-(1:10)], S$dS.Sdt[-(1:10)]))), max(c(H$dH.Hdt[-(1:10)], S$dS.Sdt[-(1:10)]))), 
         xlab="", ylab="d-1", lwd=3, cex=1, cex.lab=1)
    if(any(c(H$dH.Hdt[-(1:10)], S$dS.Sdt[-(1:10)])<0)) abline(h=0, lwd=0.5, lty=3)
    text(time[0.95*length(time)], Hgrf, labels=as.character(round(Hgrf, 3)), pos=3, xpd=T, cex=0.75)
    #title("State variable dynamics", adj=0, cex.main=2, outer = T)
    title("Specific growth rate", adj=0.05, line=0)
    lines(time, S$dS.Sdt, col="black", lwd=1)
    legend("topright", legend=c("Host", "Sym"), lwd=c(3,1), col="black", bty="n", y.intersp=1)
    
    # Symbiont to host ratio
    totSH <- S$S / H$H
    totSHf <- totSH[length(totSH)]
    plot(time, totSH, type="l", col="black", ylim=c(0, max(totSH, na.rm=T)), ylab="C-molS/C-molH", xlab="", lwd=3)
    text(time[0.95*length(time)], totSHf, labels=as.character(round(totSHf, 3)), pos=3, xpd=T, cex=0.75)
    title("Symbiont to host ratio", adj=0.05, line=0)

    # Light quenching (carbon fixation, NPQ, and excess (=ROS producing))
    plot(NA, xlim=range(time), xlab="", ylim=c(0, max(c(S$jL, S$jeL))), ylab="mol photons/C-molS/d")
    title("Light quenching", adj=0.05, line=0)
    lines(time, S$jL, col="yellow", lty=3, lwd=3) # total amount absorbed
    lines(time, S$jL - S$jeL, col="yellow", lwd=1) # amt. used in photosynthesis
    lines(time, (S$jL - S$jeL) + pmin(S$jeL, pars$jNPQ), col="yellow", lty=1, lwd=3) # amt. quenched by NPQ
    legend("topright", legend=c("Excess", "NPQ", "Photo."), lty=c(3,1,1), lwd=c(3,3,1), col="yellow", bty="n")
    
    # Symbiont biomass SU dynamics (proportions of C and N rejected from SU)
    plot(NA, xlim=range(time), xlab="", ylim=c(0, 1), ylab="Proportion of uptake")
    title("Symbiont: substrate excess", adj=0.05, line=0)
    lines(time, S$rhoC/S$jCP, col="red", lwd=3)
    lines(time, S$jNw/(H$rhoN*H$H/S$S + S$rNS), col="blue", lwd=3)
    legend("topright", legend=c("C", "N"), lwd=3, bty="n", col=c("red", "blue"))
      
    # Photosynthesis SU dynamics (proportions of light and CO2 rejected from SU)
    plot(NA, xlim=range(time), xlab="", ylim=c(0,1), ylab="Proportion of uptake")
    title("Photosynthesis: substrate excess", adj=0.05, line=0)
    lines(time, S$jeL/S$jL, col="yellow", lwd=3)
    lines(time, S$jCO2w/((H$jCO2 + H$rCH)*H$H/S$S + S$rCS), col="red", lwd=3)
    legend("topright", legend=c("Light", "DIC"), lty=c(1,1), lwd=3, bty="n", col=c("yellow", "red"))
    
    # Photosynthesis rate (carbon fixation)
    plot(NA, xlim=range(time), xlab="", ylim=c(0, max(S$jCP)), ylab="C-mol/C-molS/d")
    title("Photosynthesis rate", adj=0.05, line=0)
    lines(time, S$jCP, col="red", lwd=3)
    
    # Relative ROS production due to excitation energy in excess of carbon fixation and NPQ
    plot(NA, xlim=range(time), xlab="", ylim=c(1, max(S$cROS)), ylab="Relative to baseline")
    title("ROS production", adj=0.05, line=0)
    lines(time, S$cROS, col="orange", lwd=3)
      
    # Host biomass SU dynamics (proportions of C and N rejected from SU)
    plot(NA, xlim=range(time), xlab="", ylim=c(0, 1), ylab="Proportion of uptake")
    title("Coral: substrate excess", adj=0.05, line=0)
    lines(time, H$jCw/(S$rhoC*S$S/H$H + H$jX), col="red", lwd=3)
    lines(time, H$rhoN/(H$jN + H$jX*pars$nNX + H$rNH), col="blue", lwd=3)
    legend("topright", legend=c("C", "N"), lwd=3, bty="n", col=c("red", "blue"))
  })
}