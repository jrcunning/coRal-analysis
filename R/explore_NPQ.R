# This function plots the photosynthetic rate and relative ROS production as a function of irradiance
# for a given set of parameters (i.e., jCPm, jNPQ, kROS, k). The argument to this function should be a 
# full set of model parameters (i.e., object with structure matching output of "def_pars.R"). This script 
# is meant to show the photosynthetic performance of a given symbiont independent from a host (i.e., in 
# a state that is not CO2-limited), so the "jCO2" parameter is set at a very high value here.

sym <- function(pars, SU="pc", jL, jCO2, NPQ, NEC=NULL, kROS=NULL, photoinhibition=T) {
  if (SU=="pc") synth <- function(x, y, m) 1/((1/m)+(1/x)+(1/y)-(1/(x+y)))
  if (SU=="he") synth <- function(x, y, m) (m^-2 + x^-2 + y^-2)^(-1/2)
  with(pars, {
    # Set light range for PI curve
    time <- seq(1,100,0.1)
    jCO2 <- seq(jCO2[1], jCO2[2], along=time)
    jL <- seq(jL[1], jL[2], along=time)
    # Initialize values...
    jCP <- c(pars$jCPm, rep(NA, length(time)-1))
    jeL <- c(0, rep(NA, length(time)-1))
    jNPQ <- c(0, rep(NA, length(time)-1))
    cROS <- c(1, rep(NA, length(time)-1))
    # Time-stepping solution...
    for (t in 2:length(time)) {
      # Calculate photosynthesis rate
      if (photoinhibition==T) jCP[t] <- synth(jCO2[t], jL[t]*pars$nLC, pars$jCPm)/cROS[t-1]
      if (photoinhibition==F) jCP[t] <- synth(jCO2[t], jL[t]*pars$nLC, pars$jCPm)
      
      # Calculate light in excess of photosynthetic quenching
      jeL[t] <- max(0, jL[t] - jCP[t]/pars$nLC)
      # Calculate light energy quenched by NPQ
      if (NPQ=="none") {
        jNPQ[t] <- 0
      } else if (NPQ=="inhibited") {  #
        jNPQ[t] <- min(jeL[t], pars$kNPQ / cROS[t-1])
      } else if (NPQ=="inhibited2") {  #
        jNPQ[t] <- (pars$kNPQ^(-3)+jeL[t]^(-3))^(-1/3) / cROS[t-1]
      } else if (NPQ=="fixed") {
        jNPQ[t] <- min(jeL[t], pars$kNPQ)
      } else if (NPQ=="fixed2") {
        jNPQ[t] <- (pars$kNPQ^(-3)+jeL[t]^(-3))^(-1/3)
      } else if (NPQ=="jCP-linked") {
        jNPQ[t] <- min(jeL[t], 2 * jCP[t]/pars$nLC)
      }
      # Calculate ROS (cROS) produced due to excess light
      cROS[t] <- 1 + ((jeL[t] - jNPQ[t]) / pars$kROS)^pars$k
      if (!is.null(NEC)) {
        cROS[t] <- 1 + (max(0, jeL[t]-NEC)/pars$kROS)^pars$k
      }
      # How about:: instead of dividing jCP by cROS, divide by jeL-jNPQ? (amount of eL left after NPQ.)
    }
    # Return ROS and photosynthesis rate for plotting
    par(mfrow=c(1,1), mar=c(3,3,1,3), mgp=c(1.2,0.2,0), cex=0.7, tck=0.025, xaxs="i")
    #plot(time, jCP, xlab="Light (mol photons/C-molS/d)", ylab="Carbon fixation (mol C/C-molS/d)",
    #     main="Photosynthesis", type="l", lwd=2, col="red", ylim=c(0,3))
    #plot(time, jeL, main="jeL", type="l", lwd=2); lines(time, jeL-jNPQ, lty=2, lwd=2)
    #legend("topleft", c("after photo", "after NPQ"), lty=c(1,2), lwd=2, bty="n")
    #plot(time, jNPQ, type="l", lwd=2, main="jNPQ", ylim=c(0, max(max(jNPQ, na.rm=T), 1)))
    #plot(time, cROS, xlab="Light (mol photons/C-molS/d)", ylab="ROS (relative to baseline)",
    #     main="ROS production", type="l", lwd=2, col="orange", ylim=c(1,5))
    plot(time, jL, type="l", lwd=2, ylim=c(0, max(jL)), col="gold", ylab="molph/Cmol/d")
    lines(time, jCP/pars$nLC, col="red")
    lines(time, jCP/pars$nLC + jNPQ, lty=2, col="blue")
    par(new=T)
    plot(time, cROS, ylim=c(1,max(cROS)*1.5), type="l", col="orange", yaxt="n", ylab="")
    axis(side=4)
    mtext(side=4, "cROS (relative)", cex=0.7, line=1.5)
    legend("topleft", c("jL", "jCP", "jCP+jNPQ", "ROS"), col=c("gold", "red", "blue", "orange"), lty=c(1,1,2,1), bty="n")
  })
}




