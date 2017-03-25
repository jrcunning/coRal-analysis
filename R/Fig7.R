# Fig. 7
library(coRal)
library(dplyr)
library(scales)

time <- seq(1, 200, 0.1)  # Set time
defpars <- coRal::def_pars()  # Get default parameters
pars <- defpars

# Run 1
env1 <- coRal::init_env(time=time, L=c(25,48,2), N=c(1e-7,1e-7,0), X=c(2e-7,2e-7,0))
ss1 <- with(coRal::run_coral_ss(env=list(L=35, N=1e-7, X=2e-7), pars=pars, dt=0.1), last(S$S/H$H))
run1 <- coRal::run_coral(time=time, env=env1, pars=replace(pars, "initS", ss1))
# Run 2 - increase food
env2 <- coRal::init_env(time=time, L=c(25,48,2), N=c(1e-7,1e-7,0), X=c(1e-6,1e-6,0))
ss2 <- with(coRal::run_coral_ss(env=list(L=35, N=1e-7, X=1e-6), pars=pars, dt=0.1), last(S$S/H$H))
run2 <- coRal::run_coral(time=time, env=env2, pars=replace(pars, "initS", ss2))
# Run 3 - increase DIN
env3 <- coRal::init_env(time=time, L=c(25,48,2), N=c(4e-6,4e-6,0), X=c(2e-7,2e-7,0))
ss3 <- with(coRal::run_coral_ss(env=list(L=35, N=4e-6, X=2e-7), pars=pars, dt=0.1), last(S$S/H$H))
run3 <- coRal::run_coral(time=time, env=env3, pars=replace(pars, "initS", ss3))


# Plot
png("img/Fig7.png", width=3, height=3, units="in", res=300)
par(mfrow=c(1,1), mar=c(2,2.5,1,2.5), mgp=c(1.75,-0.2,0), tcl=0.25, xaxs="i", yaxs="i")
plot(time, env1$L, col=alpha("gold", 0.5), type="l", lwd=2, axes=F, xlab="", 
     ylab="S:H biomass", ylim=c(30,50), cex.lab=0.9)
mtext(side=1, text="Days", line=1, cex=0.9)
mtext(side=2, line=1, text=expression((C-mol~S~C-mol~H^{-1})), cex=0.5)
axis(side=4, cex.axis=0.7); mtext(side=4, text="Light", line=1, cex=0.9)
mtext(side=4, line=1.75, text=expression((mol~photons~m^{-2}~d^{-1})), cex=0.5)
par(new=T, mgp=c(2,0,0))
plot(NA, xlim=range(time), ylim=c(0.0,0.25), yaxs="i", ylab="", xlab="", cex.axis=0.7, cex.lab=0.9)
with(run1, lines(time, S$S/H$H, lty=1))
with(run2, lines(time, S$S/H$H, lty=2))
with(run3, lines(time, S$S/H$H, lty=3))
legend("bottomleft", legend=c("+DIN", "+feeding", "default"), lty=c(3,2,1), bty="n", cex=0.5)
legend("topright", legend="Light", col="gold", lty=1, bty="n", cex=0.5)
dev.off()

