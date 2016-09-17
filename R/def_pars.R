def_pars <- function() {

  # Default host parameters
  defrpars <- list(
    gammaR=0.03,  # Host specific biomass turnover rate (d-1)
    nNR=0.19,  # N:C ratio in host biomass (-)
    nNX=0.2,  # N:C ratio in prey biomass (-)
    sigmaNR=0.9,  # Proportion of host nitrogen turnover recycled (-)
    sigmaCR=0.9,  # Proportion of host carbon turnover recycled (-)
    UXm=0.1292,  # Maximum specific host feeding rate (molX/C-molR/d)
    UNm=0.048,  # Maximum specific host DIN uptake rate (molN/C-molR/d)
    #UCPm=2,  # Maximum specific host CO2 uptake/delivery rate (molCO2/C-molS/d)
    UCPpf=4.04e-3,  # Specific flux of CO2 into tissue due to passive reaciton-diffusion (Tansik et al. 2015)
    UCPaf=0.32,  # Specific delivery of CO2 to photosynthesis by coral host (Tansik et al. 2015)
    QRm=1,  # Maximum specific host growth rate (C-molR/C-molR/d)
    XKN=0.46e-6,  # Half-saturation constant for host DIN uptake (molN/L)
    XKX=20e-6,  # Half-saturation constant for host feeding (molX/L)
    #XKC=4e-4,  # Half-saturation constant for CO2 uptake (molCO2/L)
    initR=1  # Initial host biomass (C-molR)
  )  
  
  # Default symbiont parameters
  defspars <- list(
    gammaS=0.03,  # Symbiont specific biomass turnover rate (d-1),
    nNS=0.13,  # N:C ratio in symbiont biomass (-)
    nLC=0.1,  # L:C ratio in fixed carbon (=quantum yield) (mol C/mol photons)
    NPQ=40,  # no effect concentration of ROS per symbiont
    L50=40,  # amount of (eL-NPQ) that doubles ROS production
    k=1,  # exponent on ROS production
    ds=1.34,  # Symbiont specific cross-sectional area (m^2/C-molS)
    sigmaNS=0.9,  # Proportion of symbiont nitrogen turnover recylced (-)
    sigmaCS=0.9,  # Proportion of symbiont carbon turnover recycled (-)
    UCSm=2.8,  # Maximum specific photosynthate production rate (C-mol/C-molS/d)
    QSm=0.25,  # Maximum specific symbiont growth rate (C-molS/C-molS/d)
    initS=0.1  # Initial symbiont biomass (C-molS)
  )
  
  pars <- list(rpars=defrpars, spars=defspars)
  return(pars)
  
}

