# Plot environment for single model run
plot_L <- function(run) with(run, {
  plot(time, env$L, type="l", col="yellow", ylim=c(0,60), lwd=3, xlab="", ylab="mol/m2/d", main="Irradiance")
})

plot_DIN <- function(run) with(run, {
  plot(time, env$N * 10^6, type="l", col="blue", ylim=c(0,4), lwd=3, xlab="", ylab="µmol/L", main="DIN")
})

plot_gr <- function(run) with(run, {
  Sgr <- S$jSG - S$ST
  Hgr <- H$jHG - H$HT
  Hgrf <- Hgr[length(Hgr)]
  plot(time, Hgr, type="l", ylim=c(min(0, min(c(Hgr[-(1:10)], Sgr[-(1:10)]))), max(c(Hgr[-(1:10)], Sgr[-(1:10)]))),
       xlab="", ylab="d-1", lwd=3, cex=1, cex.lab=1, 
       main="Specific growth rate")
  if(any(c(Hgr[-(1:10)], Sgr[-(1:10)])<0)) abline(h=0, lwd=0.5, lty=3)
  text(time[0.95*length(time)], Hgrf, labels=as.character(round(Hgrf, 3)), pos=3, xpd=T, cex=0.75)
  lines(time, Sgr, col="black", lwd=1)
  legend("topright", legend=c("Host", "Sym"), lwd=c(3,1), col="black", bty="n", y.intersp=1)
})

plot_sh <- function(run) with(run, {
  totSH <- S$S / H$H
  totSHf <- totSH[length(totSH)]
  plot(time, totSH, type="l", col="black", ylim=c(0, max(totSH, na.rm=T)), ylab="C-molS/C-molH", xlab="", lwd=3,
       main="Symbiont:host biomass")
  text(time[0.95*length(time)], totSHf, labels=as.character(round(totSHf, 3)), pos=3, xpd=T, cex=0.75)
})

plot_Lq <- function(run) with(run, {
  plot(NA, xlim=range(time), xlab="", ylim=c(0, max(c(S$jL, S$jeL))), ylab="mol photons/C-molS/d",
       main="Light quenching")
  lines(time, S$jL, col="yellow", lty=3, lwd=3) # total amount absorbed
  lines(time, S$jL - S$jeL, col="yellow", lwd=1) # amt. used in photosynthesis
  lines(time, (S$jL - S$jeL) + pmin(S$jeL, spars$jNPQ), col="yellow", lty=1, lwd=3) # amt. quenched by NPQ
  legend("topright", legend=c("Excess", "NPQ", "Photo."), lty=c(3,1,1), lwd=c(3,3,1), col="yellow", bty="n")
})

plot_photo <- function(run) with(run, {
  plot(NA, xlim=range(time), xlab="", ylim=c(0, max(S$jCP)), ylab="C-mol/C-molS/d")
  title("Photosynthesis rate", adj=0.05, line=0)
  lines(time, S$jCP, col="red", lwd=3)
})


plot_PSU <- function(run) with(run, {
  plot(NA, xlim=range(time), xlab="", ylim=c(0,1), ylab="Proportion of uptake")
  title("Photosynthesis: substrate excess", adj=0.05, line=0)
  lines(time, S$jeL/S$jL, col="yellow", lwd=3)
  lines(time, S$eC/((H$jCO2 + H$rCH)*H$H/S$S + S$rCS), col="red", lwd=3)
  legend("topright", legend=c("Light", "DIC"), lty=c(1,1), lwd=3, bty="n", col=c("yellow", "red"))
})

plot_ROS <- function(run) with(run, {
  plot(NA, xlim=range(time), xlab="", ylim=c(1, max(S$cROS)), ylab="Relative to baseline")
  title("ROS production", adj=0.05, line=0)
  lines(time, S$cROS, col="orange", lwd=3)
})

plot_symSU <- function(run) with(run, {
  plot(NA, xlim=range(time), xlab="", ylim=c(0, 1), ylab="Proportion of uptake")
  title("Symbiont: substrate excess", adj=0.05, line=0)
  lines(time, S$rhoC/S$jCP, col="red", lwd=3)
  lines(time, S$wN/(H$rhoN*H$H/S$S + S$rNS), col="blue", lwd=3)
  legend("topright", legend=c("C", "N"), lwd=3, bty="n", col=c("red", "blue"))
})

plot_corSU <- function(run) with(run, {
  plot(NA, xlim=range(time), xlab="", ylim=c(0, 1), ylab="Proportion of uptake")
  title("Coral: substrate excess", adj=0.05, line=0)
  lines(time, H$wC/(S$rhoC*S$S/H$H + H$jX), col="red", lwd=3)
  lines(time, H$rhoN/(H$jN + H$jX*hpars$nNX + H$rNH), col="blue", lwd=3)
  legend("topright", legend=c("C", "N"), lwd=3, bty="n", col=c("red", "blue"))
})





plot_run <- function(run) {
  with(run, {
    # Set up plot
    par(mar=c(2,2,1,1))
    par(oma=c(0,0,0,0), mfrow=c(5,2), mgp=c(1.2,0.2,0), tck=0.025, lwd=1, xaxs="i", 
        cex.main=1.5, cex.axis=1, cex.lab=1)
    # External PAR
    plot(time, env$L, type="l", col="yellow", ylim=c(0,60), lwd=3, xlab="", ylab="mol/m2/d")
    #title("Environment", adj=0, cex.main=2, outer = T)
    title("Light", adj=0.05, line=0)
    # External [DIN]
    plot(time, env$N * 10^6, type="l", col="blue", ylim=c(0,4), ylab="µmol/L", xlab="", lwd=3)
    title("DIN", adj=0.05, line=0)
    
    # Plot growth rates
    Sgr <- S$jSG - S$ST
    Hgr <- H$jHG - H$HT
    Hgrf <- Hgr[length(Hgr)]
    plot(time, Hgr, type="l", ylim=c(min(0, min(c(Hgr[-(1:10)], Sgr[-(1:10)]))), max(c(Hgr[-(1:10)], Sgr[-(1:10)]))), 
         xlab="", ylab="d-1", lwd=3, cex=1, cex.lab=1)
    if(any(c(Hgr[-(1:10)], Sgr[-(1:10)])<0)) abline(h=0, lwd=0.5, lty=3)
    text(time[0.95*length(time)], Hgrf, labels=as.character(round(Hgrf, 3)), pos=3, xpd=T, cex=0.75)
    #title("State variable dynamics", adj=0, cex.main=2, outer = T)
    title("Specific growth rate", adj=0.05, line=0)
    lines(time, Sgr, col="black", lwd=1)
    legend("topright", legend=c("Host", "Sym"), lwd=c(3,1), col="black", bty="n", y.intersp=1)
    
    # Plot symbiont to host ratios
    totSH <- S$S / H$H
    totSHf <- totSH[length(totSH)]
    plot(time, totSH, type="l", col="black", ylim=c(0, max(totSH, na.rm=T)), ylab="C-molS/C-molH", xlab="", lwd=3)
    text(time[0.95*length(time)], totSHf, labels=as.character(round(totSHf, 3)), pos=3, xpd=T, cex=0.75)
    title("Symbiont to host ratio", adj=0.05, line=0)

    # UL and jeL
    plot(NA, xlim=range(time), xlab="", ylim=c(0, max(c(S$jL, S$jeL))), ylab="mol photons/C-molS/d")
    #title("Flux dynamics", adj=0, cex.main=2/0.66, outer = T)
    title("Light quenching", adj=0.05, line=0)
    lines(time, S$jL, col="yellow", lty=3, lwd=3) # total amount absorbed
    lines(time, S$jL - S$jeL, col="yellow", lwd=1) # amt. used in photosynthesis
    lines(time, (S$jL - S$jeL) + pmin(S$jeL, spars$jNPQ), col="yellow", lty=1, lwd=3) # amt. quenched by NPQ
    legend("topright", legend=c("Excess", "NPQ", "Photo."), lty=c(3,1,1), lwd=c(3,3,1), col="yellow", bty="n")
      
    # CO2 uptake per symbiont
    #jCO2 <- H$jCO2*H$H/S$S  # renormalize carbon uptake to CmolS
    #metC <- H$rCH * H$H / S$S + S$rCS
    #plot(NA, xlim=range(time), xlab="", ylim=c(0, min(3, max(c(jCO2, metC), na.rm=T))), ylab="molC/C-molS/d")
    #title("DIC uptake", adj=0.05, line=0)
    #lines(time, jCO2, col="red", lwd=3)
    #lines(time, metC, col="red", lwd=1)
    #legend("topright", legend=c("External", "Metabolic"), lty=1, lwd=c(3,1), bty="n", col="red")
    
    # Symbiont: rejection to uptake ratios
    plot(NA, xlim=range(time), xlab="", ylim=c(0, 1), ylab="Proportion of uptake")
    title("Symbiont: substrate excess", adj=0.05, line=0)
    lines(time, S$rhoC/S$jCP, col="red", lwd=3)
    lines(time, S$wN/(H$rhoN*H$H/S$S + S$rNS), col="blue", lwd=3)
    legend("topright", legend=c("C", "N"), lwd=3, bty="n", col=c("red", "blue"))
      
    # Photosynthesis: rejection to uptake ratios
    plot(NA, xlim=range(time), xlab="", ylim=c(0,1), ylab="Proportion of uptake")
    title("Photosynthesis: substrate excess", adj=0.05, line=0)
    lines(time, S$jeL/S$jL, col="yellow", lwd=3)
    lines(time, S$eC/((H$jCO2 + H$rCH)*H$H/S$S + S$rCS), col="red", lwd=3)
    legend("topright", legend=c("Light", "DIC"), lty=c(1,1), lwd=3, bty="n", col=c("yellow", "red"))
    
    # Photosynthesis rate per symbiont
    plot(NA, xlim=range(time), xlab="", ylim=c(0, max(S$jCP)), ylab="C-mol/C-molS/d")
    title("Photosynthesis rate", adj=0.05, line=0)
    lines(time, S$jCP, col="red", lwd=3)
    
    # ROS production per symbiont
    plot(NA, xlim=range(time), xlab="", ylim=c(1, max(S$cROS)), ylab="Relative to baseline")
    title("ROS production", adj=0.05, line=0)
    lines(time, S$cROS, col="orange", lwd=3)
      
    # Coral: rejection to uptake ratios
    plot(NA, xlim=range(time), xlab="", ylim=c(0, 1), ylab="Proportion of uptake")
    title("Coral: substrate excess", adj=0.05, line=0)
    lines(time, H$wC/(S$rhoC*S$S/H$H + H$jX), col="red", lwd=3)
    lines(time, H$rhoN/(H$jN + H$jX*hpars$nNX + H$rNH), col="blue", lwd=3)
    legend("topright", legend=c("C", "N"), lwd=3, bty="n", col=c("red", "blue"))
  })
}

  
