# This script initializes and runs a simulation using specified environmental inputs and model parameters.
# The arguments to this function are:
#   - time: a vector of time steps at which the model should be evaluated (units=days) (e.g., seq(0,365,0.1))
#   - env: an object with structure matching the output of "init_env.R"
#   - pars: an object with structure matching the output of "def_pars.R"


run_coral <- function(time, env, pars) {
  require(dplyr)
  # Define parallel complementary synthesizing unit function
  synth <- function(x, y, m) 1/((1/m)+(1/x)+(1/y)-(1/(x+y)))
  
  # Set initial values
  # ==================
  # Host
  H <- data_frame(
    time=time,
    H=pars$initH,
    jX=(pars$jXm * env$X[1] / (env$X[1] + pars$KX)),
    jCO2a=pars$jCO2a,
    jCO2p=pars$jCO2p,
    jCO2=jCO2a+jCO2p,
    jN=(pars$jNm * env$N[1] / (env$N[1] + pars$KN)),
    rhoN=jN,
    jCw=0,
    jHG=0.25,
    jHT=pars$jHT0,
    rNH=jHT * pars$nNH * pars$sigmaNH,
    rCH=jHT * pars$sigmaCH,
    dH.Hdt=0
  )
  H[2:nrow(H), 2:ncol(H)] <- NA
  
  # Symbiont
  S <- data_frame(
    time=time,
    S=pars$initS,
    jL=env$L[1] * pars$astar,
    jCP=max(0, synth(jL * pars$nLC, (H$jCO2p[1] + H$jCO2a[1])*H$H/S, pars$jCPm), na.rm=T),
    jeL=max(jL - jCP/pars$nLC, 0),
    jNPQ=pars$jNPQ/pars$nLC,
    jCO2w=H$jCO2*H$H/S - jCP,
    jSG=pars$jSGm,
    rhoC=jCP, 
    jNw=0,
    jST=pars$jST0,
    rNS=jST * pars$nNS * pars$sigmaNS,
    rCS=jST * pars$sigmaCS,
    cROS=1,
    dS.Sdt=1)
  S[2:nrow(S), 2:ncol(S)] <- NA
  
  # Run simulation by updating
  # ==========================
  for (t in 2:length(time)) {
    
    # Photosynthesis SU
    # =================
    # Light input flux
    S$jL[t] <- (1.256307 + 1.385969 * exp(-6.479055 * (S$S[t-1]/H$H[t-1]))) * env$L[t] * pars$astar
    # CO2 input flux
    H$jCO2p[t] <- pars$jCO2p # passive CO2 flux from the environment
    H$jCO2a[t] <- pars$jCO2a  # active delivery of CO2 from the environment (mediated by the host)
    H$jCO2[t] <- H$jCO2p[t] + H$jCO2a[t]  # Total CO2 input from the environment
    S$rCS[t] <- pars$jST0 * pars$sigmaCS  # metabolic CO2 recycled from symbiont biomass turnover
    H$rCH[t] <- H$jHT[t-1] * pars$sigmaCH  # metabolic CO2 recycled from host biomass turnover
    # Production flux (photosynthetic carbon fixation)
    S$jCP[t] <- synth(S$jL[t] * pars$nLC, (H$jCO2[t] + H$rCH[t])*H$H[t-1]/S$S[t-1] + S$rCS[t], pars$jCPm) / S$cROS[t-1]
    # Rejection flux: CO2 (wasted to the environment)
    S$jCO2w[t] <- max((H$jCO2[t] + H$rCH[t])*H$H[t-1]/S$S[t-1] + S$rCS[t] - S$jCP[t], 0)
    # Rejection flux: Light (=excitation energy not quenched by carbon fixation)
    S$jeL[t] <- max(S$jL[t] - S$jCP[t]/pars$nLC, 0)
    # Amount of excitation energy quenched by NPQ
    S$jNPQ[t] <- (pars$jNPQ * S$jCP[t])/pars$nLC
    # Scaled ROS production due to excess excitation energy (=not quenched by carbon fixation AND NPQ)
    S$cROS[t] <- 1 + (max(0, (S$jeL[t] - S$jNPQ[t])) / pars$kROS) ^ pars$k
    
    # Symbiont biomass SU
    # ===================
    # Nitrogen input flux
    S$rNS[t] <- pars$jST0 * pars$nNS * pars$sigmaNS  # Recylced N from symbiont biomass turnover.
    H$rhoN[t-1] <- H$rhoN[t-1]  # Nitrogen shared from the host (defined below, so previous time step used)
    # Carbon input flux
    S$jCP[t] <- S$jCP[t]  # Production of fixed carbon from photosynthesis SU
    # Production flux (symbiont biomass formation)
    S$jSG[t] <- synth(S$jCP[t], (H$rhoN[t-1]*H$H[t-1]/S$S[t-1] + S$rNS[t])/pars$nNS, pars$jSGm)
    # Rejection flux: carbon (surplus carbon shared with the host)
    S$rhoC[t] <- max(S$jCP[t] - S$jSG[t], 0)
    # Rejection flux: nitrogen (surplus nitrogen wasted to the environment)
    S$jNw[t] <- max(H$rhoN[t-1]*H$H[t-1]/S$S[t-1] + S$rNS[t] - pars$nNS * S$jSG[t], 0)
    
    # Symbiont biomass loss (turnover)
    S$jST[t] <- pars$jST0 * (1 + 5 * (S$cROS[t] - 1))
    
    # Host biomass SU
    # ===============
    # Food input flux (prey=both carbon and nitrogen)
    H$jX[t] <- (pars$jXm * env$X[t] / (env$X[t] + pars$KX))  # Prey uptake from the environment
    # Nitrogen input flux
    H$jN[t] <- (pars$jNm * env$N[t] / (env$N[t] + pars$KN))  # N uptake from the environment
    H$rNH[t] <- H$jHT[t-1] * pars$nNH * pars$sigmaNH  # Recycled N from host biomass turnover
    # Carbon input flux
    S$rhoC[t] <- S$rhoC[t]  # Carbon shared by the symbiont (defined above)
    # Production flux (host biomass formation)
    H$jHG[t] <- synth(S$rhoC[t]*S$S[t-1]/H$H[t-1] + H$jX[t], (H$jN[t] + pars$nNX*H$jX[t] + H$rNH[t]) / pars$nNH, pars$jHGm)
    # Rejection flux: nitrogen (surplus nitrogen shared with the symbiont)
    H$rhoN[t] <- max(H$jN[t] + pars$nNX * H$jX[t] + H$rNH[t] - pars$nNH * H$jHG[t], 0)
    # Rejection flux: carbon (surplus carbon lost to environment)
    H$jCw[t] <- max(H$jX[t] + S$rhoC[t]*S$S[t-1]/H$H[t-1] - H$jHG[t], 0)
    
    # Host biomass loss
    H$jHT[t] <- pars$jHT0
   
    # State equations
    # ===============
    # Specific growth rates (Cmol/Cmol/d)
    S$dS.Sdt[t] <- S$jSG[t] - S$jST[t]
    H$dH.Hdt[t] <- H$jHG[t] - H$jHT[t]
    # Biomass (Cmol)
    S$S[t] <- S$S[t-1] + S$dS.Sdt[t] * S$S[t-1] * (time[2] - time[1])
    H$H[t] <- H$H[t-1] + H$dH.Hdt[t] * H$H[t-1] * (time[2] - time[1])

  }
  
  # Return results
  # ==============
  return(list(time=time, env=env, pars=pars, H=H, S=S))
}

