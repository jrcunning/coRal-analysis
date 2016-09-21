# This script creates the surface plots used in Figure 2.
library(png)
source("R/def_pars.R")
source("R/init_env.R")
source("R/run_coral.R")
source("R/coral_steady_states.R")


# Define dimensions for surface plots
dim <- list(
  L=seq(from=0, to=40, length.out=21),
  N=seq(from=0, to=4e-6, length.out=21)
)

# Set up and run simulations for surface plot using default parameters
time <- seq(1,100,1)
env <- init_env(time, L=c(10,10,0), N=c(1e-7,1e-7,0), X=c(0,0,0))
pars <- def_pars()
surf <- coral_steady_states(time, env, pars, dim)

# Plot growth surface and save as png
png("img/growth.png", width=5, height=4.2, units="in", res=300)
par(mar=c(3,3,3,0), mgp=c(1.5,0.2,0), tck=0.02)
filled.contour(x=dim[[2]]*10^6, y=dim[[1]], z=t(surf$gr), nlevels=max(surf$gr) %/% 0.01,
               col=gray.colors(max(surf$gr) %/% 0.01+1, start=1, end=0.2),
               xlab="[DIN] (µmol/L)", ylab="Irradiance (mol ph/m2/d)")
title("A. Specific growth rate", adj=0)
#mtext(side=4, text="Cmol/Cmol/d")
dev.off()

# Plot growth surface and save as png
png("img/sh.png", width=5, height=4.2, units="in", res=300)
par(mar=c(3,3,3,0), mgp=c(1.5,0.2,0), tck=0.02)
filled.contour(x=dim[[2]]*10^6, y=dim[[1]], z=t(surf$sh), nlevels=max(surf$sh) %/% 0.05,
               col=gray.colors(max(surf$sh) %/% 0.05+1, start=1, end=0.2),
               xlab="[DIN] (µmol/L)", ylab="Irradiance (mol ph/m2/d)")
title("B. Symbiont:host biomass", adj=0)
dev.off()

# Plot both pngs as a two-panel figure
png1 <- readPNG("img/growth.png")
png2 <- readPNG("img/sh.png")

png("img/Fig2.png", width=7, height=2.94, units="in", res=300)
par(mfrow=c(1,2), mar=c(0,0,0,0))
plot.new()
rasterImage(png1, par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4])
plot.new()
rasterImage(png2, par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4])
dev.off()
