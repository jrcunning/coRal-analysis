# Load libraries
library(parallel)
library(doParallel)
library(foreach)

# Load functions
sapply(c("R/def_pars.R",
         "R/run_coral_ss.R"), 
       source, .GlobalEnv)

# Get default parameter values
defpars <- def_pars()

# Set up cluster for parallel processing
cl <- makeCluster(detectCores())  # Initiate cluster
registerDoParallel(cl)

# Set input values for steady state runs
input <- expand.grid(L=seq(20,55,1), initS=c(0.0001,0.1))

# Run steady states in parallel
output <- foreach(i=1:nrow(input), .combine=rbind) %dopar% {
  run <- run_coral_ss(env=list(L=input$L[i], X=2e-6, N=1e-7), pars=replace(defpars, "initS", input$initS[i]), dt=0.1)
  list(gr=last(run$H$dH.Hdt), sh=last(run$S$S/run$H$H))
}

stopCluster(cl)  # Stop cluster

# Collect results
res <- cbind(input, output)

# Create plot
png("img/newFig2.png", width=5, height=5, units="in", res=300)

plot(NA, xlim=c(20,55), ylim=c(0,0.15), xlab="Light", ylab="Steady state S:H ratio")
points(c(48,48,53,53), c(0.15,0.17,0.15,0.17), pch=c(19,1,19,1), cex=c(0.4,1,0.4,1), col=c("black","black","red","red"), xpd=T)
text(38, 0.16, labels="init. S:H", xpd=T, cex=0.7, srt=90)
text(46, c(0.15,0.17), labels=c("0.1", "0.0001"), xpd=T, pos=2, cex=0.7)
text(50, 0.20, labels=c("Steady state growth"), xpd=T, cex=0.7)
text(c(48,53), 0.185, labels=c("pos.", "neg."), xpd=T, cex=0.7)

apply(res, 1, FUN=function(x) {
  with(x, points(L, sh,
                 col=ifelse(gr>0, "black", "red"), 
                 pch=ifelse(initS==1e-01, 19, 1), 
                 cex=ifelse(initS==1e-01, 0.4, 1)))
})

dev.off()