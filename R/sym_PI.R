sym_PI <- function(pars) {
  with(pars, {
    jL <- seq(0,200,1)
    # Set carbon delivery rate (mol C / CmolS / d)
    jCO2 <- 100  # Extremely high carbon delivery - no carbon limitation
    # Calculate photosynthesis rate
    jCP <- synth(jCO2, jL*spars$nLC, spars$jCPm)
    # Calculate light in excess of photosynthetic quenching
    jeL <- pmax(0, jL - jCP/spars$nLC)
    # Calculate ROS (cROS) produced due to excess light
    cROS <- pmax(1, 1 + (pmax(0, jeL-spars$jNPQ)/(spars$kROS))^spars$k)
    # Calculate photosynthesis rate accounting for ROS-induced photoinhibition
    jCP <- jCP / cROS
    # Return ROS and photosynthesis rate for plotting
    par(mfrow=c(1,2), mar=c(2,2,1,1), mgp=c(1.2,0.2,0), cex=0.7, tck=0.025, xaxs="i")
    plot(jL, jCP, xlab="Light (mol photons/C-molS/d)", ylab="Carbon fixation (mol C/C-molS/d)",
         main="Photosynthesis", type="l", lwd=3, col="red", ylim=c(0,3))
    plot(jL, cROS, xlab="Light (mol photons/C-molS/d)", ylab="ROS (relative to baseline)",
         main="ROS production", type="l", lwd=3, col="orange", ylim=c(1,5))
  })
}

  