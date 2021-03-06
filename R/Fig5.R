# Figure 5
library(coRal)
library(dplyr)
library(scales)

# Set parameters
defpars <- coRal::def_pars()  # Get default parameters

# Set run time vector
time <- seq(1, 365, 0.1)

# Initialize environments
# Values based on Stimson et al. 1997
env1 <- coRal::init_env(time=time, L=c(20,44,2), 
                        N=c(0.14e-6,0.14e-6,0), X=c(1e-6,1e-6,0))
env2 <- coRal::init_env(time=time, L=c(20,44,2), 
                        N=c(15.14e-6,15.14e-6,0), X=c(1e-6,1e-6,0))

# Run simulations
ss1 <- with(coRal::run_coral_ss(env=list(L=32,N=0.14e-6,X=1e-6), pars=defpars, dt=0.1), dplyr::last(S$S/H$H))
run1 <- coRal::run_coral(time=time, env=env1, pars=replace(defpars, "initS", ss1))
ss2 <- with(coRal::run_coral_ss(env=list(L=32,N=15.14e-6,X=1e-6), pars=defpars, dt=0.1), dplyr::last(S$S/H$H))
run2 <- coRal::run_coral(time=time, env=env2, pars=replace(defpars, "initS", ss2))

# Create figure
png("img/Fig5.png", width=3, height=4, units="in", res=300)
layout(mat=matrix(c(1,2,2,3,3)))
par(mar=c(0,4,1,1), mgp=c(2.5,0.1,0), tcl=0.25, xaxs="i")
plot(run1$env.L ~ time, type="l", lwd=3, col=alpha("gold", 0.5), ylim=c(20,45), yaxs="i", 
     xlab="", ylab="Light", xaxt="n", bty="n", xpd=NA)
mtext(side=2, text=expression((mol~photons~m^{-2}~d^{-1})), cex=0.5, line=1.2)
par(mar=c(1.5,4,1.5,1))
plot(NA, xlim=range(time), ylim=range(0,0.3), xlab="", ylab="S:H biomass", bty="n", xaxt="n")
mtext(side=2, text=expression((C-mol~S~C-mol~H^{-1})), cex=0.5, line=1.2)
with(run1, lines(S/H ~ time))
with(run2, lines(S/H ~ time, lty=2))
par(mar=c(3,4,0,1))
plot(NA, xlim=range(time), ylim=range(0,0.15), xlab="", ylab="Specific growth", bty="n")
mtext(side=1, text="Days", line=1.5, cex=0.66)
mtext(side=2, text=expression((d^{-1})), cex=0.5, line=1.2)
with(run1, lines(dH.Hdt ~ time))
with(run2, lines(dH.Hdt ~ time, lty=2))
legend("topleft", legend=c("DIN=ambient", "DIN=enriched"), lty=c(1,2), inset=c(0.12, -0.15), xpd=NA)
dev.off()

