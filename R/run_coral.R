# This script initializes and runs a simulation using specified environmental inputs and model parameters.
# The arguments to this function are:
#   - time: a vector of time steps at which the model should be evaluated (units=days) (e.g., seq(0,365,0.1))
#   - env: an object with structure matching the output of "init_env.R"
#   - pars: an object with structure matching the output of "def_pars.R"


run_coral <- function(time, env, pars) {
  require(dplyr)
  
  # Set initial values
  # =============
  # Host
  hpars <- pars$hpars
  H <- data_frame(
    time=time,
    H=hpars$initH,
    jX=(hpars$jXm * env$X[1] / (env$X[1] + hpars$KX)),
    jCO2a=hpars$jCO2a,
    jCO2p=hpars$jCO2p,
    jCO2=jCO2a+jCO2p,
    jN=(hpars$jNm * env$N[1] / (env$N[1] + hpars$KN)),
    rhoN=jN,
    jCw=0,
    jHG=0.25,
    jHT=hpars$jHT0,
    rNH=jHT * hpars$nNH * hpars$sigmaNH,
    rCH=jHT * hpars$sigmaCH
  )
  H[2:nrow(H), 2:ncol(H)] <- NA
  
  # Symbiont
  spars <- pars$spars
  S <- data_frame(
    time=time,
    S=spars$initS,
    jL=env$L[1] * spars$astar,
    jCP=max(0, synth(jL * spars$nLC, (H$jCO2p[1] + H$jCO2a[1])*H$H/S, spars$jCPm), na.rm=T),
    jeL=max(jL - jCP/spars$nLC, 0),
    jCO2w=H$jCO2*H$H/S - jCP,
    jSG=spars$jSGm,
    rhoC=jCP, 
    jNw=0,
    jST=spars$jST0,
    rNS=jST * spars$nNS * spars$sigmaNS,
    rCS=jST * spars$sigmaCS,
    cROS=1)
  S[2:nrow(S), 2:ncol(S)] <- NA
  
  # Run simulation by updating
  # ==========================
  for (t in 2:length(time)) {
    # Symbiont biomass loss (turnover)
    S$jST[t] <- spars$jST0 * (1 + 5 * (S$cROS[t-1] - 1))
    # Host biomass loss
    H$jHT[t] <- hpars$jHT0
    # C and N recycling
    S$rNS[t] <- spars$jST0 * spars$nNS * spars$sigmaNS
    H$rNH[t] <- H$jHT[t] * hpars$nNH * hpars$sigmaNH
    S$rCS[t] <- spars$jST0 * spars$sigmaCS
    H$rCH[t] <- H$jHT[t] * hpars$sigmaCH
    # Light uptake by symbiont
    S$jL[t] <- (1.256307 + 1.385969 * exp(-6.479055 * (S$S[t-1]/H$H[t-1]))) * env$L[t] * spars$astar
    # CO2 delivery to photosynthesis SU
    H$jCO2p[t] <- hpars$jCO2p # passive CO2 flux
    H$jCO2a[t] <- hpars$jCO2a  # active delivery of CO2 by host
    H$jCO2[t] <- (H$jCO2p[t] + H$jCO2a[t])  # Total CO2 delivery
    # Photosynthetic production
    S$jCP[t] <- synth(S$jL[t] * spars$nLC, H$jCO2[t]*H$H[t-1]/S$S[t-1] + S$rCS[t] + H$rCH[t]*H$H[t-1]/S$S[t-1], spars$jCPm) / S$cROS[t-1]
    # Symbiont biomass production
    S$jSG[t] <- synth(S$jCP[t], (H$rhoN[t-1]*H$H[t-1]/S$S[t-1] + S$rNS[t])/spars$nNS, spars$jSGm)
    # Prey uptake by host
    H$jX[t] <- (hpars$jXm * env$X[t] / (env$X[t] + hpars$KX))
    # DIN uptake by host
    H$jN[t] <- (hpars$jNm * env$N[t] / (env$N[t] + hpars$KN))
    # Host biomass production
    H$jHG[t] <- synth(S$rhoC[t-1]*S$S[t-1]/H$H[t-1] + H$jX[t], (H$jN[t] + hpars$nNX*H$jX[t] + H$rNH[t]) / hpars$nNH, hpars$jHGm)
    # Symbiont sharing of surplus carbon
    S$rhoC[t] <- (S$jCP[t] - S$jSG[t])
    # Coral sharing of surplus nitrogen
    H$rhoN[t] <- H$jN[t] + hpars$nNX * H$jX[t] + H$rNH[t] - hpars$nNH * H$jHG[t]
    # Wasted nitrogen from symbiont biomass SU (lost to environment)
    S$jNw[t] <- H$rhoN[t]*H$H[t-1]/S$S[t-1] + S$rNS[t] - spars$nNS * S$jSG[t]
    # Wasted carbon from coral biomass SU (lost to environment)
    H$jCw[t] <- H$jX[t] + S$rhoC[t]*S$S[t-1]/H$H[t-1] - H$jHG[t]
    # Wasted CO2 from photosynthesis SU (lost to environment)
    S$jCO2w[t] <- H$jCO2[t]*H$H[t-1]/S$S[t-1] + H$rCH[t]*H$H[t-1]/S$S[t-1] + S$rCS[t] - S$jCP[t]
    # Excess excitation energy from photosynthesis SU (=not quenched by carbon fixation)
    S$jeL[t] <- S$jL[t] - S$jCP[t]/spars$nLC
    # Scaled ROS production due to excess excitation energy (=not quenched by carbon fixation AND NPQ)
    S$cROS[t] <- 1 + (pmax(0, (S$jeL[t] - spars$jNPQ)) / spars$kROS) ^ spars$k
    # Balance equations: change in symbiont and host biomass
    S$S[t] <- S$S[t-1] + (time[2] - time[1]) * (S$jSG[t] - S$jST[t]) * S$S[t-1]
    H$H[t] <- H$H[t-1] + (time[2] - time[1]) * (H$jHG[t] - H$jHT[t]) * H$H[t-1]
  }
  
  # Return results
  # ==============
  return(list(time=time, env=env, hpars=hpars, spars=spars, H=H, S=S))
}

