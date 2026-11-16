# Builders for test-sim.R: a model with two etas, one subject's doses and
# observations, and its fixed parameters.

sim_model <- function() {
  pmx_ode_model({
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    ka <- exp(tka)
    cp <- linCmt()
  })
}

sim_events <- function(times = c(1, 4, 12, 24)) {
  pmx_ode_events(amt = 100, cmt = "depot", time = times)
}

sim_theta <- function() c(tcl = 0, tv = 3, tka = -0.5)
