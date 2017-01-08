# Load functions
sapply(c("R/def_pars.R",
         "R/run_coral_ss.R"), 
       source, .GlobalEnv)

# Get default parameter values
defpars <- def_pars()

# Simulations with different prey N:C ratios
time <- seq(1,300,0.1)
run0.4 <- run_coral_ss(env=list(L=10, N=1e-7, X=1e-6), pars=replace(defpars, c("initS", "nNX"), c(0.0001, 0.4)), dt=0.1)
run0.2 <- run_coral_ss(env=list(L=10, N=1e-7, X=1e-6), pars=replace(defpars, "initS", 0.0001), dt=0.1)
run0.1 <- run_coral_ss(env=list(L=10, N=1e-7, X=1e-6), pars=replace(defpars, c("initS", "nNX"), c(0.0001, 0.10)), dt=0.1)
run0.05 <- run_coral_ss(env=list(L=10, N=1e-7, X=1e-6), pars=replace(defpars, c("initS", "nNX"), c(0.0001, 0.05)), dt=0.1)

# Plotting function
plotrecov <- function(run, lty) with(run, {
  lines(seq(1, nrow(S)), S$S/H$H, lty=lty)
  points(min(which(H$dH.Hdt>0)), (S$S/H$H)[min(which(H$dH.Hdt>0))], pch=19, cex=0.75)
})

# Plot
png("img/Fig8.png", width=3, height=3, units="in", res=300)
par(mfrow=c(1,1))
plot(NA, xlim=c(0,3000), ylim=c(0,1), xlab="Days", ylab="S:H ratio", xaxt="n")
axis(side=1, at=seq(0,3000,1000), labels=seq(0,3000,1000)*0.1)
plotrecov(run0.4, lty=4)
plotrecov(run0.2, lty=3)
plotrecov(run0.1, lty=2)
plotrecov(run0.05, lty=1)
legend("topleft", lty=c(4,3,2,1), legend=c("N:C=0.4", "N:C=0.2", "N:C=0.1", "N:C=0.05"), bty="n")
dev.off()
