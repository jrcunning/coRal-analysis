# Figure 8
library(coRal)
library(parallel)
library(doParallel)
library(foreach)

# Run simulations if .RData file does not exist
if (!file.exists("R/Fig8.RData")) {
  # Load default parameters
  defpars <- coRal::def_pars()
  
  # Set up cluster for parallel processing
  cl <- makeCluster(detectCores())  # Initiate cluster
  registerDoParallel(cl)
  
  # Set input values for steady state runs
  input <- expand.grid(L=seq(1,60,1), initS=c(0.0001,1))
  
  # Run steady states in parallel
  output1 <- foreach(i=1:nrow(input), .combine=rbind) %dopar% {
    run <- coRal::run_coral_ss(env=list(L=input$L[i], X=2e-7, N=1e-7), pars=replace(defpars, "initS", input$initS[i]), dt=0.1)
    list(gr=last(run$H$dH.Hdt), sh=last(run$S$S/run$H$H))
  }
  
  output2 <- foreach(i=1:nrow(input), .combine=rbind) %dopar% {
    run <- coRal::run_coral_ss(env=list(L=input$L[i], X=0, N=1e-7), pars=replace(defpars, "initS", input$initS[i]), dt=0.1)
    list(gr=last(run$H$dH.Hdt), sh=last(run$S$S/run$H$H))
  }
  
  output3 <- foreach(i=1:nrow(input), .combine=rbind) %dopar% {
    run <- coRal::run_coral_ss(env=list(L=input$L[i], X=2e-7, N=2e-6), pars=replace(defpars, "initS", input$initS[i]), dt=0.1)
    list(gr=last(run$H$dH.Hdt), sh=last(run$S$S/run$H$H))
  }
  
  output4 <- foreach(i=1:nrow(input), .combine=rbind) %dopar% {
    run <- coRal::run_coral_ss(env=list(L=input$L[i], X=4e-7, N=2e-6), pars=replace(defpars, "initS", input$initS[i]), dt=0.1)
    list(gr=last(run$H$dH.Hdt), sh=last(run$S$S/run$H$H))
  }
  
  stopCluster(cl)  # Stop cluster
  
  # Collect results
  res1 <- cbind(input, output1)
  res2 <- cbind(input, output2)
  res3 <- cbind(input, output3)
  res4 <- cbind(input, output4)
  
  save(input, res1, res2, res3, res4, file="R/Fig8.RData")
}

# Create plot
png("img/Fig8.png", width=5, height=5, units="in", res=300)
par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=c(1.5,0.5,0), tcl=-0.025)

plot(NA, xlim=range(input$L), ylim=c(0,0.5), xlab=expression("Light" ~ (mol~photons~m^{2}~d^{-1})), ylab="Steady state S:H ratio")
text(x=par("usr")[1]-par("mai")[2]*conv/1.1, y=par("usr")[4]+par("mai")[3]/2.5*conv, adj=c(0,0), 
     labels="A.", xpd=NA, font=2, cex=1)
points(c(48,48,53,53), c(0.42,0.45,0.42,0.45), pch=c(19,1,19,1), cex=c(0.4,1,0.4,1), col=c("black","black","red","red"), xpd=T)
text(38, 0.435, labels="init. S:H", xpd=T, cex=0.7, srt=90)
text(46, c(0.42,0.45), labels=c("0.1", "0.0001"), xpd=T, pos=2, cex=0.7)
text(50, 0.49, labels=c("Steady state growth"), xpd=T, cex=0.7)
text(c(48,53), 0.47, labels=c("pos.", "neg."), xpd=T, cex=0.7)
p1 <- apply(res1, 1, FUN=function(x) {
  with(x, points(L, sh,
                 col=ifelse(gr>0, "black", "red"), 
                 pch=ifelse(initS==1, 19, 1), 
                 cex=ifelse(initS==1, 0.4, 1)))
  with(x, if (sh > 0.5) {
    points(L, 0.5, pch="^", 
           col=ifelse(gr>0, "black", "red"))
    points(L, 0.49, 
           col=ifelse(gr>0, "black", "red"),
           pch=ifelse(initS==1, 19, 1), 
           cex=ifelse(initS==1, 0.4, 1))
  })
})

plot(NA, xlim=range(input$L), ylim=c(0,0.5), xlab=expression("Light" ~ (mol~photons~m^{2}~d^{-1})), ylab="Steady state S:H ratio")
text(x=par("usr")[1]-par("mai")[2]*conv/1.1, y=par("usr")[4]+par("mai")[3]/2.5*conv, adj=c(0,0), 
     labels="B.", xpd=NA, font=2, cex=1)
p2 <- apply(res2, 1, FUN=function(x) {
  with(x, points(L, sh,
                 col=ifelse(gr>0, "black", "red"), 
                 pch=ifelse(initS==1, 19, 1), 
                 cex=ifelse(initS==1, 0.4, 1)))
  with(x, if (sh > 0.5) {
    points(L, 0.5, pch="^", 
           col=ifelse(gr>0, "black", "red"))
    points(L, 0.49, 
           col=ifelse(gr>0, "black", "red"),
           pch=ifelse(initS==1, 19, 1), 
           cex=ifelse(initS==1, 0.4, 1))
  })
})

plot(NA, xlim=range(input$L), ylim=c(0,0.5), xlab=expression("Light" ~ (mol~photons~m^{2}~d^{-1})), ylab="Steady state S:H ratio")
text(x=par("usr")[1]-par("mai")[2]*conv/1.1, y=par("usr")[4]+par("mai")[3]/2.5*conv, adj=c(0,0), 
     labels="C.", xpd=NA, font=2, cex=1)
p3 <- apply(res3, 1, FUN=function(x) {
  with(x, points(L, sh,
                 col=ifelse(gr>0, "black", "red"), 
                 pch=ifelse(initS==1, 19, 1), 
                 cex=ifelse(initS==1, 0.4, 1)))
  with(x, if (sh > 0.5) {
    points(L, 0.5, pch="^", 
           col=ifelse(gr>0, "black", "red"))
    points(L, 0.49, 
           col=ifelse(gr>0, "black", "red"),
           pch=ifelse(initS==1, 19, 1), 
           cex=ifelse(initS==1, 0.4, 1))
  })
})

plot(NA, xlim=range(input$L), ylim=c(0,0.5), xlab=expression("Light" ~ (mol~photons~m^{2}~d^{-1})), ylab="Steady state S:H ratio")
text(x=par("usr")[1]-par("mai")[2]*conv/1.1, y=par("usr")[4]+par("mai")[3]/2.5*conv, adj=c(0,0), 
     labels="D.", xpd=NA, font=2, cex=1)
p4 <- apply(res4, 1, FUN=function(x) {
  with(x, points(L, sh,
                 col=ifelse(gr>0, "black", "red"), 
                 pch=ifelse(initS==1, 19, 1), 
                 cex=ifelse(initS==1, 0.4, 1)))
  with(x, if (sh > 0.5) {
    points(L, 0.5, pch="^", 
           col=ifelse(gr>0, "black", "red"))
    points(L, 0.49, 
           col=ifelse(gr>0, "black", "red"),
           pch=ifelse(initS==1, 19, 1), 
           cex=ifelse(initS==1, 0.4, 1))
  })
})
dev.off()

