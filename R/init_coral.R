init_coral <- function(pars) {
  require(dplyr)
  
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
  
  return(list(R=R, S=S))
}
