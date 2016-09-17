run_coral <- function(time, env, pars) {
  
  require(dplyr)
  
  # Get initial values
  # Host initialization
  rpars <- pars$rpars
  R <- data_frame(
    time=time,
    R=rpars$initR,
    UX=(rpars$UXm * env$X[1] / (env$X[1] + rpars$XKX)),
    UCPa=rpars$UCPaf,
    UCPp=rpars$UCPpf,
    UCPt=UCPa+UCPp,
    UN=(rpars$UNm * env$N[1] / (env$N[1] + rpars$XKN)),
    rhoN=UN,
    wC=0,
    QR=0.25,
    TR=rpars$gammaR,
    rNR=TR * rpars$nNR * rpars$sigmaNR,
    rCR=TR * rpars$sigmaCR
  )
  R[2:nrow(R), 2:ncol(R)] <- NA
  
  # Symbiont initialization
  spars <- pars$spars
  S <- data_frame(
    time=time,
    S=spars$initS,
    UL=env$L[1] * spars$ds,
    UCS=max(0, synth(UL * spars$nLC, (R$UCPp[1] + R$UCPa[1])*R$R/S, spars$UCSm), na.rm=T),
    eL=max(UL - UCS/spars$nLC, 0),
    eC=R$UCPt*R$R/S - UCS,
    QS=spars$QSm,
    rhoC=UCS, 
    wN=0,
    TS=spars$gammaS,
    rNS=TS * spars$nNS * spars$sigmaNS,
    rCS=TS * spars$sigmaCS,
    dp=1)
  S[2:nrow(S), 2:ncol(S)] <- NA
  
  # Run simulation by updating
  for (t in 2:length(time)) {
    # Symbiont biomass loss
    # Make turnover rate increase twice the dp (ROS) flux? 5 times!!!!
    S$TS[t] <- spars$gammaS + spars$gammaS * 5 * (S$dp[t-1] - 1) #+ 5 * (S$dp[t-1] - 1) # REMOVE TIMES 5 -- should be 5 * (dp - 1)
    # Coral biomass loss
    R$TR[t] <- rpars$gammaR
    # Biomass recycling
    S$rNS[t] <- S$TS[t] * spars$nNS * spars$sigmaNS
    R$rNR[t] <- R$TR[t] * rpars$nNR * rpars$sigmaNR
    S$rCS[t] <- S$TS[t] * spars$sigmaCS
    R$rCR[t] <- R$TR[t] * rpars$sigmaCR
    # Light uptake
    #S$UL[t] <- env$L[t] * (1 + 2 * exp(-3.8 * (S$S[t-1]/R$R[t-1]))) * spars$ds   #- old
    #S$UL[t] <- (0.6775 + 3.6947 * exp(-43.8762 * (S$S[t-1]/R$R[t-1]))) * env$L[t] * spars$ds  # aboral data only
    #S$UL[t] <- (0.9977 + 2.4820 * exp(-39.0402 * (S$S[t-1]/R$R[t-1]))) * env$L[t] * spars$ds  # aboral and oral data
    #S$UL[t] <- (0.9815 + (3 - 0.9815) * exp(-33.2613 * max(0, (S$S[t-1]/R$R[t-1])))) * env$L[t] * spars$ds  # aboral and oral with y-intercept=3
    #S$UL[t] <- 1.667 * exp(-3.739 * (S$S[t-1]/R$R[t-1])) * env$L[t] * spars$ds  # nl2 -- no intercept
    S$UL[t] <- (1.256307 + 1.385969 * exp(-6.479055 * (S$S[t-1]/R$R[t-1]))) * env$L[t] * spars$ds  # from Marcelino data
    # DIC flux to photosynthesis SU
    R$UCPp[t] <- rpars$UCPpf # passive CO2 flux due to reaction-diffusion
    R$UCPa[t] <- rpars$UCPaf  # active delivery of CO2 by host to symbiont photosynthesis SU
    R$UCPt[t] <- (R$UCPp[t] + R$UCPa[t])
    # Photosynthetic production
    S$UCS[t] <- synth(S$UL[t] * spars$nLC, R$UCPt[t]*R$R[t-1]/S$S[t-1] + S$rCS[t] + R$rCR[t]*R$R[t-1]/S$S[t-1], spars$UCSm) / S$dp[t-1]
    # Symbiont biomass production
    S$QS[t] <- synth(S$UCS[t], (R$rhoN[t-1]*R$R[t-1]/S$S[t-1] + S$rNS[t])/spars$nNS, spars$QSm)
    # Heterotrophic food uptake
    R$UX[t] <- (rpars$UXm * env$X[t] / (env$X[t] + rpars$XKX))
    # DIN uptake
    R$UN[t] <- (rpars$UNm * env$N[t] / (env$N[t] + rpars$XKN))
    # Coral biomass production
    R$QR[t] <- synth(S$rhoC[t-1]*S$S[t-1]/R$R[t-1] + R$UX[t], (R$UN[t] + rpars$nNX*R$UX[t] + R$rNR[t]) / rpars$nNR, rpars$QRm)
    # Symbiont sharing of surplus carbon
    S$rhoC[t] <- (S$UCS[t] - S$QS[t])
    # Coral sharing of surplus nitrogen
    R$rhoN[t] <- R$UN[t] + rpars$nNX * R$UX[t] + R$rNR[t] - rpars$nNR * R$QR[t]
    # Wasted nitrogen from symbiont
    S$wN[t] <- R$rhoN[t]*R$R[t-1]/S$S[t-1] + S$rNS[t] - spars$nNS * S$QS[t]
    # Wasted carbon from coral
    R$wC[t] <- R$UX[t] + S$rhoC[t]*S$S[t-1]/R$R[t-1] - R$QR[t]
    # Excess DIC from photosynthesis
    S$eC[t] <- R$UCPt[t]*R$R[t-1]/S$S[t-1] + R$rCR[t]*R$R[t-1]/S$S[t-1] + S$rCS[t] - S$UCS[t]
    # Excess light energy from photosynthesis
    S$eL[t] <- S$UL[t] - S$UCS[t]/spars$nLC
    # Damage potential (ROS) due to excess light energy
    S$dp[t] <- 1 + (pmax(0, (S$eL[t] - spars$NPQ)) / spars$L50) ^ spars$k
    # Balance equations: change in biomass
    S$S[t] <- S$S[t-1] + (time[2] - time[1]) * (S$QS[t] - S$TS[t]) * S$S[t-1]
    R$R[t] <- R$R[t-1] + (time[2] - time[1]) * (R$QR[t] - R$TR[t]) * R$R[t-1]
  }
  
  # Return results
  return(list(time=time, env=env, rpars=rpars, spars=spars, R=R, S=S))
}

