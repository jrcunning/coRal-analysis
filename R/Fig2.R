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
input <- expand.grid(L=seq(20,50,1), initS=c(0.0001,1))

# Run steady states in parallel
output <- foreach(i=1:nrow(input), .combine=rbind) %dopar% {
  run <- run_coral_ss(env=list(L=input$L[i], X=1e-7, N=1e-7), pars=replace(defpars, "initS", input$initS[i]), dt=0.1)
  list(gr=last(run$H$dH.Hdt), sh=last(run$S$S/run$H$H))
}

stopCluster(cl)  # Stop cluster

# Collect results
res <- cbind(input, output)

# Create plot
png("img/Fig2.png", width=5, height=5, units="in", res=300)

plot(NA, xlim=c(20,50), ylim=c(0,0.20), xlab="Light", ylab="Steady state S:H ratio")
points(c(43,43,48,48), c(0.17,0.18,0.17,0.18), pch=c(19,1,19,1), cex=c(0.4,1,0.4,1), col=c("black","black","red","red"), xpd=T)
text(41, 0.20, labels="init. S:H", xpd=T, cex=0.7, pos=2)
text(41, c(0.17,0.18), labels=c("1", "0.0001"), xpd=T, pos=2, cex=0.7)
text(46, 0.20, labels=c("Steady state growth"), xpd=T, cex=0.7)
text(c(43,48), 0.19, labels=c(">0", "<0"), xpd=T, cex=0.7)

apply(res, 1, FUN=function(x) {
  with(x, points(L, sh,
                 col=ifelse(gr>0, "black", "red"), 
                 pch=ifelse(initS==1, 19, 1), 
                 cex=ifelse(initS==1, 0.4, 1)))
  with(x, if (sh > 0.2) {
    points(L, 0.205, pch="^", 
           col=ifelse(gr>0, "black", "red"))
    points(L, 0.20, 
           col=ifelse(gr>0, "black", "red"),
           pch=ifelse(initS==1, 19, 1), 
           cex=ifelse(initS==1, 0.4, 1))
  })
})

dev.off()
