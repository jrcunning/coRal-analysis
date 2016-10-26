# This script initializes and runs a simulation using specified environmental inputs and model parameters.
# The arguments to this function are:
#   - time: a vector of time steps at which the model should be evaluated (units=days) (e.g., seq(0,365,0.1))
#   - env: an object with structure matching the output of "init_env.R"
#   - pars: an object with structure matching the output of "def_pars.R"


run_coral_ss <- function(env, pars) {
  require(dplyr)
  # Define parallel complementary synthesizing unit function
  synth <- function(x, y, m) 1/((1/m)+(1/x)+(1/y)-(1/(x+y)))
  
  # Initialize environment
  
  # Set initial values
  # ==================
  # Host
  H <- data_frame(
    H=pars$initH,
    jX=(pars$jXm * env$X / (env$X + pars$KX)),
    jN=(pars$jNm * env$N / (env$N + pars$KN)),
    rhoN=jN,
    jeC=10,
    jCO2=pars$kCO2 * jeC,
    jHG=0.25,
    jHT=pars$jHT0,
    rNH=jHT * pars$nNH * pars$sigmaNH,
    rCH=jHT * pars$sigmaCH,
    dH.Hdt=0
  )
  
  # Symbiont
  S <- data_frame(
    S=pars$initS,
    jL=env$L * pars$astar,
    jCP=max(0, synth(jL * pars$yCL, H$jCO2[1]*H$H/S, pars$jCPm), na.rm=T),
    jeL=max(jL - jCP/pars$yCL, 0),
    jNPQ=pars$kNPQ,
    jCO2w=H$jCO2*H$H/S - jCP,
    jSG=pars$jSGm/10,
    rhoC=jCP, 
    jNw=0,
    jST=pars$jST0,
    rNS=jST * pars$nNS * pars$sigmaNS,
    rCS=jST * pars$sigmaCS,
    cROS=1,
    dS.Sdt=0)
  
  
  # Run simulation by updating
  # ==========================
  grss <- FALSE
  shss <- FALSE
  t <- 2
  
  while (grss==FALSE | shss==FALSE) {
    
    # Photosynthesis
    # ==============
    # Light input flux
    S[t,]$jL <- (1.256307 + 1.385969 * exp(-6.479055 * (S$S[t-1]/H$H[t-1]))) * env$L * pars$astar
    # CO2 input flux
    S[t,]$rCS <- pars$jST0 * pars$sigmaCS  # metabolic CO2 recycled from symbiont biomass turnover
    #H[t,]$rCH <- H$jHT[t-1] * pars$sigmaCH  # metabolic CO2 recycled from host biomass turnover
    H[t,]$rCH <- H$jHT[t-1] * S$S[t-1]/H$H[t-1]  #
    H[t,]$jCO2 <- pars$kCO2 * H$jeC[t-1]  # carbon not used in host biomass is used to activate CCM's that deliver CO2 to photosynthesis
    # Production flux (photosynthetic carbon fixation)
    S[t,]$jCP <- synth(S[t,]$jL * pars$yCL, (H[t,]$jCO2 + H[t,]$rCH)*H$H[t-1]/S$S[t-1] + S[t,]$rCS, pars$jCPm) / S$cROS[t-1]
    # Rejection flux: CO2 (wasted to the environment)
    S[t,]$jCO2w <- max((H[t,]$jCO2 + H[t,]$rCH)*H$H[t-1]/S$S[t-1] + S[t,]$rCS - S[t,]$jCP, 0)
    # Rejection flux: excess light energy not quenched by carbon fixation
    S[t,]$jeL <- max(S[t,]$jL - S[t,]$jCP/pars$yCL, 0)
    # Amount of excess light energy quenched by NPQ
    S[t,]$jNPQ <- (pars$kNPQ^(-1)+S[t,]$jeL^(-1))^(-1/1)
    # Scaled ROS production due to excess excitation energy (=not quenched by carbon fixation AND NPQ)
    S[t,]$cROS <- 1 + ((S[t,]$jeL - S[t,]$jNPQ) / pars$kROS)^pars$k
    
    # Symbiont biomass
    # ================
    # Nitrogen input flux
    S[t,]$rNS <- pars$jST0 * pars$nNS * pars$sigmaNS  # Recylced N from symbiont biomass turnover.
    H$rhoN[t-1] <- H$rhoN[t-1]  # Nitrogen shared from the host (defined below, so previous time step used)
    # Carbon input flux
    S[t,]$jCP <- S[t,]$jCP  # Production of fixed carbon from photosynthesis SU
    # Production flux (symbiont biomass formation)
    S[t,]$jSG <- synth(S[t,]$jCP, (H$rhoN[t-1]*H$H[t-1]/S$S[t-1] + S[t,]$rNS)/pars$nNS, pars$jSGm)
    # Rejection flux: carbon (surplus carbon shared with the host)
    S[t,]$rhoC <- max(S[t,]$jCP - S[t,]$jSG, 0)
    # Rejection flux: nitrogen (surplus nitrogen wasted to the environment)
    S[t,]$jNw <- max(H$rhoN[t-1]*H$H[t-1]/S$S[t-1] + S[t,]$rNS - pars$nNS * S[t,]$jSG, 0)
    # Symbiont biomass loss (turnover)
    S[t,]$jST <- pars$jST0 * (1 + pars$b * (S[t,]$cROS - 1))
    
    # Host biomass
    # ============
    # Food input flux (prey=both carbon and nitrogen)
    H[t,]$jX <- (pars$jXm * env$X / (env$X + pars$KX))  # Prey uptake from the environment
    # Nitrogen input flux
    H[t,]$jN <- (pars$jNm * env$N / (env$N + pars$KN))  # N uptake from the environment
    H[t,]$rNH <- H$jHT[t-1] * pars$nNH * pars$sigmaNH  # Recycled N from host biomass turnover
    # Carbon input flux
    S[t,]$rhoC <- S[t,]$rhoC  # Carbon shared by the symbiont (defined above)
    # Production flux (host biomass formation)
    H[t,]$jHG <- synth(S[t,]$rhoC*S$S[t-1]/H$H[t-1] + H[t,]$jX, (H[t,]$jN + pars$nNX*H[t,]$jX + H[t,]$rNH) / pars$nNH, pars$jHGm)
    # Rejection flux: nitrogen (surplus nitrogen shared with the symbiont)
    H[t,]$rhoN <- max(H[t,]$jN + pars$nNX * H[t,]$jX + H[t,]$rNH - pars$nNH * H[t,]$jHG, 0)
    # Rejection flux: carbon -- given back to symbiont as CO2 input to photosynthesis
    H[t,]$jeC <- max(H[t,]$jX + S[t,]$rhoC*S$S[t-1]/H$H[t-1] - H[t,]$jHG, 0)
    # Host biomass loss
    H[t,]$jHT <- pars$jHT0
    
    # State equations
    # ===============
    # Specific growth rates (Cmol/Cmol/d)
    S[t,]$dS.Sdt <- S[t,]$jSG - S[t,]$jST
    H[t,]$dH.Hdt <- H[t,]$jHG - H[t,]$jHT
    # Biomass (Cmol)
    S[t,]$S <- S$S[t-1] + S[t,]$dS.Sdt * S$S[t-1] * 1
    H[t,]$H <- H$H[t-1] + H[t,]$dH.Hdt * H$H[t-1] * 1
   
    # Test if steady state has been reached
    if (t > 21) {
      grss <- ifelse(abs(H$dH.Hdt[t] - H$dH.Hdt[t-20]) < 0.00001, T, F)
      shss <- ifelse(abs(S$S[t]/H$H[t] - S$S[t-20]/H$H[t-20]) < 0.00001, T, F)
    }
    # Test if system is oscillating but with negative growth (e.g., L=1.2, N=2.6e-6, X=0)
    if (t > 1000) {
      if (any(H$dH.Hdt[(t-500):t] < 0)) {
        grss <- T; shss <- T
        H$dH.Hdt[t] <- 0
      }
    } 

    # Increment time
    t <- t + 1
  }
  
  # Return results
  # ==============
  return(list(env=env, pars=pars, H=H, S=S))
}

