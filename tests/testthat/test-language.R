# The language (design §6.1): every construct compiles, and what it computes
# is what R computes for the same statements on the same values. R is the
# oracle because the language is R syntax; rxode2 is the second one for
# the promise that what compiles here compiles there.


test_that("each expression construct computes R's value", {
  inputs <- list(t = 3.5, y = c(x = 4.25), par = c(a = 0.7, b = 1.9))
  for (e in construct_cases()) {
    text <- paste0("y <- ", e, "\nd/dt(x) <- -y * x")
    expect_matches_r(text, inputs, e)
  }
})

test_that("statements compute R's values: =, if/else, specials, cmt()", {
  inputs <- list(t = 2, y = c(depot = 30, central = 5, peri = 1),
                 par = c(a = 0.7, b = 1.9, c = -0.4))
  models <- list(
    eq = "k = a * 2\nd/dt(central) = -k * central",
    if_else = "if (a < b) {\n k <- a\n} else {\n k <- b\n}\nd/dt(central) <- -k * central",
    if_only = "k <- 1\nif (c > 0) k <- 2\nd/dt(central) <- -k * central",
    if_state = "if (central > 3) {\n k <- a\n} else {\n k <- b\n}\nd/dt(central) <- -k * central",
    nested = "if (a > 0) {\n if (b > 2) k <- 1 else k <- 2\n} else {\n k <- 3\n}\nd/dt(central) <- -k * central",
    deriv_in_if = "if (a < b) {\n d/dt(central) <- -a * central\n} else {\n d/dt(central) <- -b * central\n}",
    specials = paste("cmt(depot)\ncmt(central)",
                     "alag(depot) <- a\nf(depot) <- b / 2\nrate(central) <- 3 * b",
                     "dur(central) <- a + 1\ncentral(0) <- 10 * b",
                     "d/dt(depot) <- -a * depot\nd/dt(central) <- a * depot - b * central",
                     sep = "\n"),
    synonyms = "lag(depot) <- a\nF(depot) <- b\nd/dt(depot) <- -a * depot",
    reassign = "k <- a\nk <- k * 2\nk <- k + b\nd/dt(central) <- -k * central",
    three = "d/dt(depot) <- -a * depot\nd/dt(central) <- a * depot - b * central + c * peri\nd/dt(peri) <- -c * peri + b * central"
  )
  for (k in names(models)) expect_matches_r(models[[k]], inputs, k)
})

test_that("every corpus model without linCmt() computes R's values", {
  skip_under_gctorture()
  texts <- corpus_model_texts()
  ran <- 0L
  for (name in names(texts)) {
    text <- texts[[name]]
    if (uses_lincmt(text)) next
    m <- tryCatch(pmx_ode_model(text), pmx_ode_error_compile = function(e) NULL)
    if (is.null(m)) next
    for (seed in 1:3) {
      inputs <- model_inputs(m, seed)
      # gamma, lgamma and ^ may differ from R's in the last place.
      expect_matches_r(text, inputs, paste(name, seed), tolerance = 1e-12)
    }
    ran <- ran + 1L
  }
  expect_gte(ran, 30L)
})

test_that("linCmt() becomes rxode2's compartments and their equations", {
  par <- c(tcl = log(2), tv = log(20), tka = log(0.8), tq = log(4),
           tvp = log(60), tq2 = log(1), tvp2 = log(100))
  y <- c(depot = 50, central = 30, peripheral1 = 10, peripheral2 = 4)
  cl <- 2; v <- 20; ka <- 0.8; q <- 4; vp <- 60; q2 <- 1; vp2 <- 100
  m <- pmx_ode_model({
    cl <- exp(tcl); v <- exp(tv); ka <- exp(tka); q <- exp(tq)
    vp <- exp(tvp); q2 <- exp(tq2); vp2 <- exp(tvp2)
    cp <- linCmt()
  })
  expect_identical(m$states, c("depot", "central", "peripheral1", "peripheral2"))
  expect_identical(m$class, "closed-form")
  got <- run_blocks(m, 0, y, par[m$params])
  k <- cl / v; k12 <- q / v; k21 <- q / vp; k13 <- q2 / v; k31 <- q2 / vp2
  expect_equal(unname(got$rhs), c(
    -ka * y[["depot"]],
    ka * y[["depot"]] - k * y[["central"]] - k12 * y[["central"]] +
      k21 * y[["peripheral1"]] - k13 * y[["central"]] + k31 * y[["peripheral2"]],
    k12 * y[["central"]] - k21 * y[["peripheral1"]],
    k13 * y[["central"]] - k31 * y[["peripheral2"]]
  ), tolerance = 1e-14)
  expect_equal(got$post[["cp"]], y[["central"]] / v, tolerance = 1e-15)

  # Micro-constants, one compartment, no depot.
  m <- pmx_ode_model("v <- exp(tv)\nk <- exp(tk)\ncp <- linCmt()")
  expect_identical(m$states, "central")
  got <- run_blocks(m, 0, c(central = 7), c(tv = log(20), tk = log(0.1)))
  expect_equal(got$rhs[["central"]], -0.1 * 7, tolerance = 1e-15)
  expect_identical(m$lincmt$ncmt, 1L)
  expect_true(m$lincmt$micro)
})

test_that("states, parameters and outputs are named in order of appearance", {
  m <- pmx_ode_model({
    ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl / v * central
    cp <- central / v
  })
  expect_identical(m$states, c("depot", "central"))
  expect_identical(m$params, c("tka", "tcl", "eta.cl", "tv"))
  expect_identical(m$captures, c("ka", "cl", "v", "cp"))
  expect_identical(m$hoisted, c("ka", "cl", "v"))
  expect_snapshot(print(m), transform = function(x) sub("(native|interp) tier", "<tier>", x))
})

test_that("a model can be given as braces, a function or text", {
  text <- "d/dt(central) <- -k * central\ncp <- central"
  a <- pmx_ode_model({
    d/dt(central) <- -k * central
    cp <- central
  })
  b <- pmx_ode_model(function() {
    d/dt(central) <- -k * central
    cp <- central
  })
  c <- pmx_ode_model(text)
  expect_identical(a$ir, c$ir)
  expect_identical(b$ir, c$ir)
  held <- text
  expect_identical(pmx_ode_model(held)$ir, c$ir)
})

test_that("what zupmxode compiles, rxode2 compiles too", {
  skip_on_cran()
  skip_under_gctorture()
  skip_if_not_installed("rxode2")
  texts <- c(
    unlist(corpus_model_texts()),
    vapply(construct_cases(), function(e) paste0("y <- ", e, "\nd/dt(x) <- -y * x"), "")
  )
  for (text in texts) {
    ok <- tryCatch({
      pmx_ode_model(text)
      TRUE
    }, pmx_ode_error_compile = function(e) FALSE)
    if (!ok) next
    # rxModelVars() parses and translates the model as rxode2() does,
    # without the C compiler that would take 89 of the suite's seconds.
    compiled <- tryCatch({
      suppressMessages(rxode2::rxModelVars(text))
      TRUE
    }, error = function(e) FALSE)
    expect_true(compiled, label = paste("rxode2 accepts", text))
  }
})
