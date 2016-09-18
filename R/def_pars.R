def_pars <- function() {

  # Default host parameters
  defhpars <- list(
    jHT0=0.03,  # Host specific biomass turnover rate (d-1)
    nNH=0.19,  # N:C ratio in host biomass (-)
    nNX=0.2,  # N:C ratio in prey biomass (-)
    sigmaNH=0.9,  # Proportion of host nitrogen turnover recycled (-)
    sigmaCH=0.9,  # Proportion of host carbon turnover recycled (-)
    jXm=0.1292,  # Maximum specific host feeding rate (molX/C-molH/d)
    jNm=0.048,  # Maximum specific host DIN uptake rate (molN/C-molH/d)
    jCO2p=4.04e-3,  # Specific flux of CO2 into tissue due to passive reaciton-diffusion (Tansik et al. 2015)
    jCO2a=0.32,  # Specific delivery of CO2 to photosynthesis by coral host (Tansik et al. 2015)
    jHGm=1,  # Maximum specific host growth rate (C-molH/C-molH/d)
    KN=0.46e-6,  # Half-saturation constant for host DIN uptake (molN/L)
    KX=20e-6,  # Half-saturation constant for host feeding (molX/L)
    #XKC=4e-4,  # Half-saturation constant for CO2 uptake (molCO2/L)
    initH=1  # Initial host biomass (C-molH)
  )  
  
  # Default symbiont parameters
  defspars <- list(
    jST0=0.03,  # Symbiont specific biomass turnover rate (d-1),
    nNS=0.13,  # N:C ratio in symbiont biomass (-)
    nLC=0.1,  # L:C ratio in fixed carbon (=quantum yield) (mol C/mol photons)
    jNPQ=40,  # no effect concentration of ROS per symbiont
    kROS=40,  # amount of (jeL-jNPQ) that doubles ROS production
    k=1,  # exponent on ROS production
    astar=1.34,  # Symbiont specific cross-sectional area (m^2/C-molS)
    sigmaNS=0.9,  # Proportion of symbiont nitrogen turnover recylced (-)
    sigmaCS=0.9,  # Proportion of symbiont carbon turnover recycled (-)
    jCPm=2.8,  # Maximum specific photosynthate production rate (C-mol/C-molS/d)
    jSGm=0.25,  # Maximum specific symbiont growth rate (C-molS/C-molS/d)
    initS=0.1  # Initial symbiont biomass (C-molS)
  )
  
  pars <- list(hpars=defhpars, spars=defspars)
  return(pars)
  
}

