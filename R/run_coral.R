run_coral <- function(time, env, pars) {
  require(dplyr)
  
  # Get initial values 
  # Host initialization
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
    wC=0,
    jHG=0.25,
    HT=hpars$jHT0,
    rNH=HT * hpars$nNH * hpars$sigmaNH,
    rCH=HT * hpars$sigmaCH
  )
  H[2:nrow(H), 2:ncol(H)] <- NA
  
  # Symbiont initialization
  spars <- pars$spars
  S <- data_frame(
    time=time,
    S=spars$initS,
    jL=env$L[1] * spars$astar,
    jCP=max(0, synth(jL * spars$nLC, (H$jCO2p[1] + H$jCO2a[1])*H$H/S, spars$jCPm), na.rm=T),
    jeL=max(jL - jCP/spars$nLC, 0),
    eC=H$jCO2*H$H/S - jCP,
    jSG=spars$jSGm,
    rhoC=jCP, 
    wN=0,
    ST=spars$jST0,
    rNS=ST * spars$nNS * spars$sigmaNS,
    rCS=ST * spars$sigmaCS,
    cROS=1)
  S[2:nrow(S), 2:ncol(S)] <- NA
  
  # Run simulation by updating
  for (t in 2:length(time)) {
    # Symbiont biomass loss
    # Make turnover rate increase twice the cROS (ROS) flux? 5 times!!!!
    S$ST[t] <- spars$jST0 * (1 + 5 * (S$cROS[t-1] - 1)) #+ 5 * (S$cROS[t-1] - 1) # REMOVE TIMES 5 -- should be 5 * (cROS - 1)
    # Coral biomass loss
    H$HT[t] <- hpars$jHT0
    # Biomass recycling
    S$rNS[t] <- spars$jST0 * spars$nNS * spars$sigmaNS
    H$rNH[t] <- H$HT[t] * hpars$nNH * hpars$sigmaNH
    S$rCS[t] <- spars$jST0 * spars$sigmaCS
    H$rCH[t] <- H$HT[t] * hpars$sigmaCH
    # Light uptake
    #S$jL[t] <- env$L[t] * (1 + 2 * exp(-3.8 * (S$S[t-1]/H$H[t-1]))) * spars$astar   #- old
    #S$jL[t] <- (0.6775 + 3.6947 * exp(-43.8762 * (S$S[t-1]/H$H[t-1]))) * env$L[t] * spars$astar  # aboral data only
    #S$jL[t] <- (0.9977 + 2.4820 * exp(-39.0402 * (S$S[t-1]/H$H[t-1]))) * env$L[t] * spars$astar  # aboral and oral data
    #S$jL[t] <- (0.9815 + (3 - 0.9815) * exp(-33.2613 * max(0, (S$S[t-1]/H$H[t-1])))) * env$L[t] * spars$astar  # aboral and oral with y-intercept=3
    #S$jL[t] <- 1.667 * exp(-3.739 * (S$S[t-1]/H$H[t-1])) * env$L[t] * spars$astar  # nl2 -- no intercept
    S$jL[t] <- (1.256307 + 1.385969 * exp(-6.479055 * (S$S[t-1]/H$H[t-1]))) * env$L[t] * spars$astar  # from Marcelino data
    # DIC flux to photosynthesis SU
    H$jCO2p[t] <- hpars$jCO2p # passive CO2 flux due to reaction-diffusion
    H$jCO2a[t] <- hpars$jCO2a  # active delivery of CO2 by host to symbiont photosynthesis SU
    H$jCO2[t] <- (H$jCO2p[t] + H$jCO2a[t])
    # Photosynthetic production
    S$jCP[t] <- synth(S$jL[t] * spars$nLC, H$jCO2[t]*H$H[t-1]/S$S[t-1] + S$rCS[t] + H$rCH[t]*H$H[t-1]/S$S[t-1], spars$jCPm) / S$cROS[t-1]
    # Symbiont biomass production
    S$jSG[t] <- synth(S$jCP[t], (H$rhoN[t-1]*H$H[t-1]/S$S[t-1] + S$rNS[t])/spars$nNS, spars$jSGm)
    # Heterotrophic food uptake
    H$jX[t] <- (hpars$jXm * env$X[t] / (env$X[t] + hpars$KX))
    # DIN uptake
    H$jN[t] <- (hpars$jNm * env$N[t] / (env$N[t] + hpars$KN))
    # Coral biomass production
    H$jHG[t] <- synth(S$rhoC[t-1]*S$S[t-1]/H$H[t-1] + H$jX[t], (H$jN[t] + hpars$nNX*H$jX[t] + H$rNH[t]) / hpars$nNH, hpars$jHGm)
    # Symbiont sharing of surplus carbon
    S$rhoC[t] <- (S$jCP[t] - S$jSG[t])
    # Coral sharing of surplus nitrogen
    H$rhoN[t] <- H$jN[t] + hpars$nNX * H$jX[t] + H$rNH[t] - hpars$nNH * H$jHG[t]
    # Wasted nitrogen from symbiont
    S$wN[t] <- H$rhoN[t]*H$H[t-1]/S$S[t-1] + S$rNS[t] - spars$nNS * S$jSG[t]
    # Wasted carbon from coral
    H$wC[t] <- H$jX[t] + S$rhoC[t]*S$S[t-1]/H$H[t-1] - H$jHG[t]
    # Excess DIC from photosynthesis
    S$eC[t] <- H$jCO2[t]*H$H[t-1]/S$S[t-1] + H$rCH[t]*H$H[t-1]/S$S[t-1] + S$rCS[t] - S$jCP[t]
    # Excess light energy from photosynthesis
    S$jeL[t] <- S$jL[t] - S$jCP[t]/spars$nLC
    # Damage potential (ROS) due to excess light energy
    S$cROS[t] <- 1 + (pmax(0, (S$jeL[t] - spars$jNPQ)) / spars$kROS) ^ spars$k
    # Balance equations: change in biomass
    S$S[t] <- S$S[t-1] + (time[2] - time[1]) * (S$jSG[t] - S$ST[t]) * S$S[t-1]
    H$H[t] <- H$H[t-1] + (time[2] - time[1]) * (H$jHG[t] - H$HT[t]) * H$H[t-1]
  }
  
  # Return results
  return(list(time=time, env=env, hpars=hpars, spars=spars, H=H, S=S))
}

