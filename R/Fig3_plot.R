library(reshape2)
library(dplyr)

source("R/def_pars.R")
source("R/init_env.R")
source("R/run_coral.R")
source("R/plot_results.R")

# Set run time vector
time <- seq(1, 365, 0.5)  # if single model run, use time input

# Initialize environment
env <- init_env(time=time, 
                L=c(5,30,2),  # Varying light
                N=c(1e-7,1e-7,0),  # Low external DIN
                X=c(0,0,0))  # No feeding

# Set parameters
pars <- def_pars()  # Get default parameters 

# Run simulation
run <- run_coral(time=time, env=env, pars=pars)

# Make figure
png("img/Fig3.png", height=6, width=3, units="in", res=300)
plot_run2(run)
dev.off()
