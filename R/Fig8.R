# Load functions
sapply(c("R/def_pars.R",
         "R/init_env.R",
         "R/run_coral.R"), 
       source, .GlobalEnv)

# Set parameters and time
defpars <- def_pars()  # Get default parameters
time <- seq(0, 250, 0.1)

# set up limitation calculation
limcoef <- function(x,y,m) log(pmin(x, m)/pmin(y, m))
#limcoef <- function(x,y,m) log(x/y)

# Run no overshoot
env1 <- init_env(time=time, L=c(20,20,0), N=c(1e-7,1e-7,0), X=c(1e-7,1e-7,0))
run1 <- run_coral(time=time, env=env1, pars=replace(defpars, c("initS", "nNX"), c(0.001, 0.1)))
run1 <- within(run1, {
  pl <- limcoef((H$jCO2 + H$rCH)*H$H/S$S + S$rCS, S$jL * pars$yCL, pars$jCPm)
  sl <- limcoef(S$jCP, (H$rhoN*H$H/S$S + S$rNS)/pars$nNS, pars$jSGm)
  hl <- limcoef(S$rhoC*S$S/H$H + H$jX, (H$jN + pars$nNX*H$jX + H$rNH) / pars$nNH, pars$jHGm)
})
# Run overshoot
env2 <- init_env(time=time, L=c(20,20,0), N=c(1e-7,1e-7,0), X=c(1e-7,1e-7,0))
run2 <- run_coral(time=time, env=env2, pars=replace(defpars, c("initS", "nNX"), c(0.001, 0.175)))
run2 <- within(run2, {
  pl <- limcoef((H$jCO2 + H$rCH)*H$H/S$S + S$rCS, S$jL * pars$yCL, pars$jCPm)
  sl <- limcoef(S$jCP, (H$rhoN*H$H/S$S + S$rNS)/pars$nNS, pars$jSGm)
  hl <- limcoef(S$rhoC*S$S/H$H + H$jX, (H$jN + pars$nNX*H$jX + H$rNH) / pars$nNH, pars$jHGm)
})
# Run larger overshoot
env3 <- init_env(time=time, L=c(20,20,0), N=c(1e-7,1e-7,0), X=c(1e-7,1e-7,0))
run3 <- run_coral(time=time, env=env3, pars=replace(defpars, c("initS", "nNX"), c(0.001, 0.22)))
run3 <- within(run3, {
  pl <- limcoef((H$jCO2 + H$rCH)*H$H/S$S + S$rCS, S$jL * pars$yCL, pars$jCPm)
  sl <- limcoef(S$jCP, (H$rhoN*H$H/S$S + S$rNS)/pars$nNS, pars$jSGm)
  hl <- limcoef(S$rhoC*S$S/H$H + H$jX, (H$jN + pars$nNX*H$jX + H$rNH) / pars$nNH, pars$jHGm)
})

# Set initial growth to 0 for plotting purposes
run1$H$dH.Hdt[1] <- 0
run2$H$dH.Hdt[1] <- 0
run3$H$dH.Hdt[1] <- 0

# Plot
png("img/Fig8.png", width=6, height=3.2, units="in", res=300)
layout(mat=matrix(c(1,2,3,4,5,6), ncol=3))
par(mgp=c(1.2,0.2,0), cex=0.66, tcl=0.25, xaxs="i", yaxs="i")
with(run1, {
  par(mar=c(1.25,2.5,2,0.5))
  plot(S$S/H$H ~ time, type="l", lty=3, ylim=c(0,0.6), bty="n", xaxt="n", xlab="", ylab="S:H biomass")
  legend("topright", bty="n", lty=c(0,3,1), lwd=1.5, col=c(NA, "black","black"), cex=c(0.8,0.6,0.6),
         legend=c("Host growth", "negative", "positive"), y.intersp=c(0,1,1), seg.len=1)
  lines(time[H$dH.Hdt>0], (S$S/H$H)[H$dH.Hdt>0])
  par(mar=c(0,0.5,2.5,0)); title("A.", adj=0)
  par(mar=c(2.5,2.5,0.75,0.5))
  plot(time, pl, type="l", ylim=c(-4,4), xlim=range(time), col="gray", bty="n", xlab="Days", ylab="Limitation coefficient")
  lines(time, sl, col="green")
  lines(time, hl, col="brown")
  abline(h=0, lty=2)
})
legend("topleft", bty="n", xpd=NA, inset=c(0.1, -0.15), seg.len=1,
       legend=c(">0 = N-limitation", "Symbiont", "Host", "<0 = C-limitation"), cex=c(0.6,0.8,0.8,0.6), pt.cex=1.25, 
       lty=c(0,1,1,0), lwd=1.5, col=c("black", "green", "brown", "black"), pch=c(NA,NA,NA,NA), y.intersp=c(0,0.5,0.7,1.1))
legend("topleft", bty="n", xpd=NA, inset=c(0.5,-0.125), seg.len=1,
       legend=c(">0 = Light-limitation", "Photosynthesis", "<0 = CO2-limitation"), cex=c(0.6,0.8,0.6), pt.cex=1.25, 
       lty=c(0,1,0), lwd=1.5, col=c("black", "gray", "black"), pch=c(NA,NA,NA), y.intersp=c(0,0.5,1))
with(run2, {
  par(mar=c(1.25,2.5,2,0.5))
  plot(S$S/H$H ~ time, type="l", lty=3, ylim=c(0,0.6), bty="n", xaxt="n", xlab="", ylab="S:H biomass")
  lines(time[H$dH.Hdt>0], (S$S/H$H)[H$dH.Hdt>0])
  par(mar=c(0,0.5,2.5,0)); title("B.", adj=0)
  par(mar=c(2.5,2.5,0.75,0.5))
  plot(time, pl, type="l", ylim=c(-4,4), xlim=range(time), col="gray", bty="n", xlab="Days", ylab="Limitation coefficient")
  lines(time, sl, col="green")
  lines(time, hl, col="brown")
  abline(h=0, lty=2)
})
with(run3, {
  par(mar=c(1.25,2.5,2,0.5))
  plot(S$S/H$H ~ time, type="l", lty=3, ylim=c(0,0.6), bty="n", xaxt="n", xlab="", ylab="S:H biomass")
  lines(time[H$dH.Hdt>0], (S$S/H$H)[H$dH.Hdt>0])
  par(mar=c(0,0.5,2.5,0)); title("C.", adj=0)
  par(mar=c(2.5,2.5,0.75,0.5))
  plot(time, pl, type="l", ylim=c(-4,4), xlim=range(time), col="gray", bty="n", xlab="Days", ylab="Limitation coefficient")
  lines(time, sl, col="green")
  lines(time, hl, col="brown")
  abline(h=0, lty=2)
})
dev.off()

