time <- seq(1, 200, 0.1)  # Set time
defpars <- def_pars()  # Get default parameters
pars <- defpars

# Run 1
env1 <- init_env(time=time, L=c(20,39,2), N=c(1e-7,1e-7,0), X=c(3e-6,3e-6,0))
ss1 <- with(run_coral_ss(env=list(L=29.6, N=1e-7, X=3e-6), pars=pars), last(S$S/H$H))
run1 <- run_coral(time=time, env=env1, pars=replace(pars, "initS", ss1))
# Run 2 - increase food
env2 <- init_env(time=time, L=c(20,39,2), N=c(1e-7,1e-7,0), X=c(6e-6,6e-6,0))
ss2 <- with(run_coral_ss(env=list(L=29.6, N=1e-7, X=6e-6), pars=pars), last(S$S/H$H))
run2 <- run_coral(time=time, env=env2, pars=replace(pars, "initS", ss2))
# Run 3 - increase DIN
env3 <- init_env(time=time, L=c(20,39,2), N=c(1e-6,1e-6,0), X=c(3e-6,3e-6,0))
ss3 <- with(run_coral_ss(env=list(L=29.6, N=1e-6, X=3e-6), pars=pars), last(S$S/H$H))
run3 <- run_coral(time=time, env=env3, pars=replace(pars, "initS", ss3))


# Plot
png("img/Fig6.png", width=3, height=3, units="in", res=300)
par(mfrow=c(1,1), mar=c(2,2,1,2), mgp=c(1,0,0), tcl=0.25, xaxs="i", yaxs="i")
plot(time, env1$L, col=alpha("gold", 0.5), type="l", lwd=2, axes=F, xlab="Days", ylab="S:H biomass ratio", ylim=c(25,40), cex.lab=0.9)
axis(side=4, cex.axis=0.75); mtext(side=4, text="Light (mol ph/m2/d)", line=1, cex=0.9)
par(new=T)
plot(NA, xlim=range(time), ylim=c(0.0,0.16), yaxs="i", ylab="", xlab="", cex.axis=0.75, cex.lab=0.9)
with(run1, lines(time, S$S/H$H, lty=1))
with(run2, lines(time, S$S/H$H, lty=2))
with(run3, lines(time, S$S/H$H, lty=3))
legend("bottomleft", legend=c("+DIN", "+feeding", "default"), lty=c(3,2,1), bty="n", cex=0.5)
legend("topright", legend="Light", col="gold", lty=1, bty="n", cex=0.5)
dev.off()

