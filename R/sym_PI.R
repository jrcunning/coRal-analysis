# This function plots the photosynthetic rate and relative ROS production as a function of irradiance
# for a given set of parameters (i.e., jCPm, jNPQ, kROS, k). The argument to this function should be a 
# full set of model parameters (i.e., object with structure matching output of "def_pars.R"). This script 
# is meant to show the photosynthetic performance of a given symbiont independent from a host (i.e., in 
# a state that is not CO2-limited), so the "jCO2" parameter is set at a very high value here.

sym_PI <- function(pars) {
  synth <- function(x, y, m) 1/((1/m)+(1/x)+(1/y)-(1/(x+y)))
  with(pars, {
    # Set light range for PI curve
    time <- seq(1,100,0.01)
    jL <- seq(0,250,along=time)
    # Set carbon delivery rate (mol C / CmolS / d)
    jCO2 <- 1000  # Extremely high carbon delivery - no carbon limitation
    # Initialize values...
    jCP <- c(pars$jCPm, rep(NA, length(time)-1))
    jeL <- c(0, rep(NA, length(time)-1))
    jNPQ <- c(0, rep(NA, length(time)-1))
    cROS <- c(1, rep(NA, length(time)-1))
    # Time-stepping solution...
    for (t in 2:length(time)) {
      # Calculate photosynthesis rate
      jCP[t] <- synth(jCO2, jL[t]*pars$nLC, pars$jCPm)/cROS[t-1]
      # Calculate light in excess of photosynthetic quenching
      jeL[t] <- max(0, jL[t] - jCP[t]/pars$nLC)
      # Calculate light energy quenched by NPQ
      jNPQ[t] <- (pars$kNPQ^(-3)+jeL[t]^(-3))^(-1/3)
      # Calculate ROS (cROS) produced due to excess light
      cROS[t] <- 1 + ((jeL[t] - jNPQ[t]) / pars$kROS)^pars$k
    }
    # Return ROS and photosynthesis rate for plotting
    par(mfrow=c(2,2), mar=c(2,2,1,1), mgp=c(1.2,0.2,0), cex=0.7, tck=0.025, xaxs="i")
    plot(jL, jCP, xlab="Light (mol photons/C-molS/d)", ylab="Carbon fixation (mol C/C-molS/d)",
         main="Photosynthesis", type="l", lwd=3, col="red", ylim=c(0,3))
    plot(jL, jeL, main="jeL")
    plot(jL, jNPQ, main="jNPQ")
    plot(jL, cROS, xlab="Light (mol photons/C-molS/d)", ylab="ROS (relative to baseline)",
         main="ROS production", type="l", lwd=3, col="orange", ylim=c(1,5))
  })
}


  