# This script creates the surface plots used in Figure 2.
library(foreach)
library(parallel)
library(doParallel)
library(reshape2)

source("R/def_pars.R")
source("R/init_env.R")
source("R/run_coral.R")
source("R/coral_steady_states.R")
source("R/coral_steady_states_par.R")


# Set time vector and default parameters
time <- seq(1,120,1)
pars <- def_pars()
#pars$jCO2a <- 1  # having higher carbon delivery from host reduces S:H ratios and causes photosynthesis to be light-limited more frequently.

# Set values of L and N at which to get steady state values
at <- expand.grid(L=seq(from=0, to=40, length.out=21),
                  N=seq(from=0, to=2e-6, length.out=21))

# Set up cluster for parallel processing
cl <- makeCluster(detectCores())  # Initiate cluster
registerDoParallel(cl)

ptm <- proc.time()
# Run simulations in parallel for each combination of L and N values
steady_states <- foreach(i=1:nrow(at), .combine='rbind', .packages='dplyr') %dopar% {
  L <- at[i,]$L; N <- at[i,]$N
  env <- init_env(time=time, L=c(L,L,0), N=c(N,N,0), X=c(0,0,0))
  run <- run_coral(time=time, env=env, pars=pars)
  ss <- lapply(run[c("H", "S")], function(x) x[length(time), ])
  gr <- ss$H$dH.Hdt
  sh <- ss$S$S/ss$H$H
  hl <- with(ss, log(S$rhoC*S$S/H$H / (H$jN+H$rNH)/pars$nNH))
  sl <- with(ss, log((H$jCO2+H$rCH)*H$H/S$S / (H$rhoN*H$H/S$S + S$rNS)))
  ee <- with(ss, max(0, S$jL - (S$jCP/run$pars$nLC + run$pars$jNPQ)))
  pl <- with(ss, log(((H$jCO2+H$rCH)*H$H/S$S + S$rCS)/(S$jL*run$pars$nLC)))
  data.frame(gr=gr, sh=sh, hl=hl, sl=sl, ee=ee, pl=pl)
}
ss <- cbind(at, steady_states)
save(ss, file = "data/ss.Rdata")

stopCluster(cl)  # Stop cluster
ptm - proc.time()  # 2601 simulations with 100 time steps each = 427 sec
                   # 10201 simulations with 100 time steps each = 1670 sec

# For runs that did NOT reach a steady state of positive growth, set all results to zero
ss[which(ss$gr < 0), c("gr", "sh", "hl", "sl", "ee", "pl")] <- NA
# For runs where steady state value of sh > 1, set sh to 1
ss[which(ss$sh > 1), "sh"] <- 1


# Reshape steady state results into matrices for plotting cont
sh <- acast(ss, L~N, value.var="sh")
gr <- acast(ss, L~N, value.var="gr")
hl <- acast(ss, L~N, value.var="hl")
sl <- acast(ss, L~N, value.var="sl")
ee <- acast(ss, L~N, value.var="ee")
pl <- acast(ss, L~N, value.var="pl")

# Build plotting function
imagef <- function(r, bin, main, col) {
  par(mar=c(3,3,3,6), mgp=c(1,0.5,0))
  if (any(r<0, na.rm=T) & any(r>0, na.rm=T)) {
    nb <- rev(seq(0, -1*min(r[which(r!=0)], na.rm=T)+bin, bin))*-1  # negative portion of breaks
    pb <- seq(0, max(r[which(r!=0)], na.rm=T)+bin, bin)
    breaks <- c(nb, pb[-1])
  } else {
    breaks <- seq(min(r[which(r!=0)], na.rm=T), max(r, na.rm=T)+bin, bin)
  }
  ncolors <- length(breaks) - 1
  if (col=="grayscale") {
    colors <- gray.colors(n=ncolors, start=0.9, end=0.1)
  } else {
    dir <- which.max(table(breaks>0))
    nn <- max(table(breaks>0))*2
    rc1 <- colorRampPalette(colors=c("red", "gray90"), space="Lab")(nn/2)
    rc1 <- rc1[-(length(rc1))]
    rc2 <- colorRampPalette(colors=c("gray90", "gray30"), space="Lab")((nn/2)+1)[-1]
    colors <- c(rc1, rc2)
    if (dir==1) {
      colors <- colors[1:ncolors]
    } else {
      colors <- colors[(length(colors)-ncolors+1):length(colors)]
    }
  }
  image(t(r), breaks=breaks, col=colors, xaxt="n", yaxt="n", main=main, adj=0)
  mtext(side=1, text="[DIN] (µmol/L)", line=2)
  axis(side=1, at=seq(0,1,0.1), 
       labels=quantile(seq(as.numeric(first(colnames(r))),
                           as.numeric(last(colnames(r))),
                           length.out=100),
                       probs=seq(0,1,0.1)))
  mtext(side=2, text="Irradiance (mol ph/m2/d)", line=2)
  axis(side=2, at=seq(0,1,0.1), 
       labels=quantile(seq(as.numeric(first(rownames(r))),
                           as.numeric(last(rownames(r))),
                           length.out=100),
                       probs=seq(0,1,0.1)))
  par(new=T, mar=c(3,2,1,0))
  
  # Plot color key
  w <- 0.05
  xl <- par("usr")[2] * (1 + w)
  yb <- par("usr")[3]
  h <- (par("usr")[4] - par("usr")[3])/ncolors
  for(i in 1:ncolors) {
    rect(xl, yb+h*(i-1), xl+w, yb+h*i, col=colors[i], xpd=NA)
  }
  axis(side=4, at=seq(par("usr")[3], par("usr")[4], h), labels=round(as.numeric(as.character(breaks)), 2), 
       pos=xl+w, padj=0.5, hadj=0, las=1, tck=0)
}

par(mfrow=c(2,2))
imagef(gr, bin=0.01, main="A. Specific growth", col="grayscale")
imagef(sh, bin=0.1, main="B. Symbiont:host biomass", col="grayscale")
imagef(hl, bin=0.25, main="C. Host N-limitation", col="grayscale")
#imagef(ee, bin=10, main="D. ")
imagef(pl, bin=0.2, main="E. Light- vs. CO2-limitation", col="RdGy")






