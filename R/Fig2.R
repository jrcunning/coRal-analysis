# Load functions
sapply(c("R/def_pars.R",
         "R/run_steady_states.R", 
         "R/plot_steady_states.R"), 
       source, .GlobalEnv)

# Set parameters
defpars <- def_pars()  # Get default parameters

# Set values of L and N at which to get steady state values
at <- expand.grid(L=seq(from=0, to=60, length.out=41),
                  N=seq(from=0, to=4e-6, length.out=41))

# Run to steady state and plot results
ssdat <- run_steady_states(pars=replace(defpars, "initS", 1), at=at, food=0, runtime=T)
plot_steady_states(ssdat, png="img/Fig2.png")
