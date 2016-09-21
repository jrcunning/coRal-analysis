coral_steady_states <- function(time, env, pars, at) {

  # Create empty matrices to populate with steady state values (gr=growth, sh=S/H ratio)
  gr <- matrix(NA, length(dim[[1]]), length(dim[[2]]), dimnames=dim)
  sh <- matrix(NA, length(dim[[1]]), length(dim[[2]]), dimnames=dim)
  
  # Calculate steady states for each combination of points
  for (x in 1:length(dim[[1]])) {
    # Create constant vector for dimension 1 at each evaluation point
    dimx <- seq(dim[[1]][x], dim[[1]][x], along=time)
    env1 <- replace(env, match(names(dim)[1], names(env)), list(dimx))
    for (y in 1:length(dim[[2]])) {
      # Create constant vector for dimension 2 at each evaluation point
      dimy <- seq(dim[[2]][y], dim[[2]][y], along=time)
      env2 <- replace(env1, match(names(dim)[2], names(env1)), list(dimy))
      # Run model to steady state
      res <- run_coral(time, env2, pars)
      # Calculate final growth rates and S/H ratios
      gr[x,y] <- res$H$dH.Hdt[length(res$H$dH.Hdt)]
      sh[x,y] <- (res$S$S/res$H$H)[length(res$S$S)]
    }
  }
  
  # Replace any negative growth rates with zero for both responses
  sh[gr<0] <- 0
  gr[gr<0] <- 0
  # Set any NA values to zero
  sh[is.na(sh)] <- 0
  gr[is.na(gr)] <- 0
  
  # Return results
  return(list(time=time, env=env, pars=pars, dim=dim, gr=gr, sh=sh))
  
}