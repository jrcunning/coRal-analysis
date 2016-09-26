# Load functions
sapply(c("R/def_pars.R",
         "R/run_steady_states.R", 
         "R/plot_steady_states.R"), 
       source, .GlobalEnv)

# Set up simulations
time <- seq(1,100,1)  # Set time vector
pars <- def_pars()  # Get default parameters
# Set values of L and N at which to get steady state values
at <- expand.grid(L=seq(from=0, to=40, length.out=41),
                  N=seq(from=0, to=4e-6, length.out=41))

# Get steady state values and plot results to png
ssdat <- run_steady_states(time=time, pars=pars, at=at)
plot_steady_states(ssdat, png="img/Fig2.png")
