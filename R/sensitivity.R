library(parallel)
library(doParallel)
library(foreach)
source("R/def_pars.R")
source("R/run_coral.R")
source("R/init_env.R")

# Sensitivity analysis

# Create function to measure sensitivity of responses to specified change in single parameter
sens <- function(env, pars, par, change) {
  # Run 1 
  run1 <- run_coral_ss(env=env, pars=pars)
  run1ss <- lapply(run1[c("H", "S")], function(x) x[nrow(run1$H), ])
  run1gr <- run1ss$H$dH.Hdt
  run1sh <- run1ss$S$S/run1ss$H$H
  # Run 2 - with changed parameter
  run2 <- run_coral_ss(env=env, pars=replace(pars, par, with(pars, get(par))*change))
  run2ss <- lapply(run2[c("H", "S")], function(x) x[nrow(run2$H), ])
  run2gr <- run2ss$H$dH.Hdt
  run2sh <- run2ss$S$S/run2ss$H$H
  # Calculate relative change in steady state response
  grchange <- run2gr / run1gr
  shchange <- run2sh / run1sh
  cROSchange <- run2ss$S$cROS / run1ss$S$cROS
  # Return relative changes
  return(list(grchange=grchange, shchange=shchange, cROSchange=cROSchange))
}


# PLOT FUNCTION
sens.plot <- function(response, pars, cols, lines="response") {
  # Collect data for given parameters
  dat <- list()
  for (par in pars) {
    dat[[par]] <- get(paste(par, ".sens", sep=""))
  }
  # Create plot
  plot(NA, xlim=range(dat[[1]]["change"]), ylim=range(dat[[1]]["change"]),
       xlab="Proportional change in parameter", ylab="Proportional change in response")
  # Add lines (and points if model moves into negative growth stable state (=NA in *.sens))
  for (par in pars) {
    for (e in levels(dat[[par]]$env)) {
      with(subset(dat[[par]], env==e), {
        for (r in response) {
          lty <- ifelse(lines=="response", match(r, levels(factor(response))), match(e, levels(dat[[par]]$env)))
          lines(change, get(r), lty=lty, col=cols[match(par, pars)])
          minl <- ifelse(any(is.na(grchange[change < 1])), min(change[!is.na(grchange) & change > min(dat[[1]]["change"])]), NA)
          maxl <- ifelse(any(is.na(grchange[change > 1])), max(change[!is.na(grchange) & change < max(dat[[1]]["change"])]), NA)
          thresh <- na.omit(c(minl, maxl))
          points(thresh, get(r)[change %in% thresh], col=cols[match(par, pars)])
        }
      })
    }
  }
 
  # Add legends
  legend("topright", legend=pars, lty=1, col=cols, bty="n", cex=0.7)
  if (lines=="response") {
    legend("topleft", legend=c("Growth", "S:H biomass"), lty=1:nlevels(factor(response)), bty="n", cex=0.7)
  } else {
    legend("topleft", legend=levels(dat[[1]]$env), lty=1:nlevels(dat[[1]]$env), bty="n", cex=0.7)
  }
}
