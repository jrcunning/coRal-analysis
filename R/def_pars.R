# This function returns the default values of all model parameters.

def_pars <- function() {
  return(list(
    jHT0=0.03,  # Host specific biomass turnover rate (d^-1)
    nNH=0.19,  # N:C ratio in host biomass (-)
    nNX=0.2,  # N:C ratio in prey biomass (-)
    sigmaNH=0.9,  # Proportion of host nitrogen turnover recycled (-)
    sigmaCH=0.9,  # Proportion of host carbon turnover recycled (-)
    jXm=0.1292,  # Maximum specific host feeding rate (molX/CmolH/d)
    jNm=0.048,  # Maximum specific host DIN uptake rate (molN/CmolH/d)
    jHGm=1,  # Maximum specific host growth rate (CmolH/CmolH/d)
    kCO2=5,  # Rate of host CCM's (molCO2/molC/d)
    KN=1e-6,  # Half-saturation constant for host DIN uptake (molN/L)
    KX=20e-6,  # Half-saturation constant for host feeding (molX/L)
    initH=1,  # Initial host biomass (CmolH)
    jST0=0.03,  # Symbiont specific biomass turnover rate (d^-1)
    nNS=0.13,  # N:C ratio in symbiont biomass (-)
    yCL=0.1,  # L:C ratio in fixed carbon (=quantum yield) (molC/mol ph)
    kNPQ=60,  # capacity of non-photochemical quenching (mol ph/CmolS/d)
    # calculated as twice the max photosynthesis rate, converted from C to photons...
    # A value of 40 was calculated from Gorbunov et al. 2001
    kROS=60,  # amount of excess light beyond NPQ capacity (e.g., jeL-jNPQ) that doubles ROS production relative to baseline (mol ph/CmolS/d)
    k=1,  # exponent on ROS production (-)
    astar=1.34,  # Symbiont specific cross-sectional area (m^2/C-molS)
    sigmaNS=0.9,  # Proportion of symbiont nitrogen turnover recylced (-)
    sigmaCS=0.9,  # Proportion of symbiont carbon turnover recycled (-)
    jCPm=2.8,  # Maximum specific photosynthate production rate (Cmol/CmolS/d)
    jSGm=0.25,  # Maximum specific symbiont growth rate (CmolS/CmolS/d)
    initS=0.1,  # Initial symbiont biomass (CmolS)
    b=5  # Scaling parameter for bleaching response
  ))
}

