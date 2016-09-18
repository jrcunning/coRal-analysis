# This function returns a list of vectors that define the external irradiance (L), 
# DIN concentration (N), and prey concentration (X) to be used as inputs for a model simulation.
# The arguments to this function are:
#   - time: a vector of time steps at which the model should be evaluated (units=days) (e.g., seq(0,365,0.1))
#   - L, N, and X (light, DIN, and prey, respectively):
#      - each of these arguments should be a vector with length=3, where the first element is the
#        minimum value, the second element is the maximum value, and the third element is the functional
#        form, where 1=linear increase, 2=linear decrease, and 3=sinusoid. The sinusoid function has a period 
#        of one year.

init_env <- function(time, L, N, X) {
  require(scales)
  # Light
  L <- if (L[3]==0) {
    seq(L[1], L[2], along=time)
  } else if (L[3]==1) {
    seq(L[2], L[1], along=time)
  } else if (L[3]==2) {
    0.5 * (L[2] - L[1]) * sin(0.0172*time) + L[1] + 0.5 * (L[2] - L[1])
  } else if (L[3]==3) {
    f <- 0.5 * (L[2] - L[1]) * sin(0.0172*time)
    f <- rescale(f, to=c(L[1]/0.4167, L[2]))
    ff <- rescale(1 + (0.01 / (1 + exp(0.05*(time-250)))), to=c(0.4167,1))
    f * ff
  }
  
  # DIN
  N <- if (N[3]==0) {
    seq(N[1], N[2], along=time)
  } else if (N[3]==1) {
    seq(N[2], N[1], along=time)
  } else {
    0.5 * (N[2] - N[1]) * sin(0.0172*time) + N[1] + 0.5 * (N[2] - N[1])
  }
  # Prey
  X <- if (X[3]==0) {
    seq(X[1], X[2], along=time)
  } else if (X[3]==1) {
    seq(X[2], X[1], along=time)
  } else {
    0.5 * (X[2] - X[1]) * sin(0.0172*time) + X[1]*10^-6 + 0.5 * (X[2] - X[1])
  }
  # Set environment specifications
  env <- list(L=L, N=N, X=X)
  return(env)
}


