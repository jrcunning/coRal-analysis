sym_PI <- function(pars) {
  with(pars, {
    UL <- seq(0,200,1)
    # Set carbon delivery rate (mol C / CmolS / d)
    UCP <- 100  # Extremely high carbon delivery - no carbon limitation
    # Calculate photosynthesis rate
    UCS <- synth(UCP, UL*spars$nLC, spars$UCSm)
    # Calculate light in excess of photosynthetic quenching
    eL <- pmax(0, UL - UCS/spars$nLC)
    # Calculate ROS (dp) produced due to excess light
    dp <- pmax(1, 1 + (pmax(0, eL-spars$NPQ)/(spars$L50))^spars$k)
    # Calculate photosynthesis rate accounting for ROS-induced photoinhibition
    UCS <- UCS / dp
    # Return ROS and photosynthesis rate for plotting
    par(mfrow=c(1,2), mar=c(2,2,1,1), mgp=c(1.2,0.2,0), cex=0.7, tck=0.025, xaxs="i")
    plot(UL, UCS, xlab="Light (mol photons/C-molS/d)", ylab="Carbon fixation (mol C/C-molS/d)",
         main="Photosynthesis", type="l", lwd=3, col="red", ylim=c(0,3))
    plot(UL, dp, xlab="Light (mol photons/C-molS/d)", ylab="ROS (relative to baseline)",
         main="ROS production", type="l", lwd=3, col="orange", ylim=c(1,5))
  })
}

  