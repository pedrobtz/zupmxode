# an eta named like an event column is refused

    Code
      pmx_ode_sim(m, sim_events(), theta = sim_theta(), omega = c(ii = 0.1))
    Condition
      Error in `pmx_ode_sim()`:
      ! `omega` names `ii`, which an event table reads as its own column.

