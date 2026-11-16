# A solver failure (design §8): NA for the subject from where it failed,
# a warning naming the subject and the time, the other subjects intact.

test_that("a subject whose model breaks down is NA from there; the others are solved", {
  # The derivative is NaN once t passes tend: no solver can go on. (A
  # finite-time blow-up such as dx/dt = x^2 is not a test both solvers
  # see: LSODA can step over the pole and report success.)
  m <- pmx_ode_model("d/dt(x) <- -x * sqrt(tend - t)")
  ev <- rbind(cbind(pmx_ode_events(amt = 1, cmt = "x", time = c(0.5, 2, 3), id = 1), tend = 1),
              cbind(pmx_ode_events(amt = 1, cmt = "x", time = c(0.5, 2, 3), id = 2), tend = 10))
  for (method in names(solver_codes)) {
    w <- NULL
    got <- withCallingHandlers(
      pmx_ode_solve(with_solver(m, method = method), ev, double()),
      pmx_ode_warning_solver = function(cnd) {
        w <<- cnd
        invokeRestart("muffleWarning")
      }
    )
    expect_false(is.null(w), label = method)
    expect_identical(w$id, 1)
    # Dormand-Prince stops where its steps fail, just past 1; LSODA runs on
    # with NaN to the end of the interval it was given: (0.5, 2] either way.
    expect_true(w$time >= 0.5 && w$time <= 2, label = paste(method, w$time))
    # x(t) = exp(-2/3 (tend^1.5 - (tend - t)^1.5)) while t <= tend.
    exact <- function(t, tend) exp(-2 / 3 * (tend^1.5 - (tend - t)^1.5))
    one <- got[got$id == 1, ]
    expect_equal(one$x[1], exact(0.5, 1), tolerance = 1e-5)
    expect_true(all(is.na(one$x[2:3])), label = method)
    two <- got[got$id == 2, ]
    expect_equal(two$x, exact(c(0.5, 2, 3), 10), tolerance = 1e-5)
  }
})
