# A saved model keeps only its IR (D18), so readRDS() is how bytes the
# package did not just make come back in: a restored IR is verified before
# anything is lowered from it.

test_that("a saved and restored model gives the same doubles", {
  hb <- hb_one_cpt()
  for (tier in available_tiers()) {
    m <- new_model(hb$ir, tier)
    case <- hb$cases[[1]]
    before <- hb$run(hb_eval_model(m), case$par, case$y)
    r <- save_restore(m)
    expect_identical(r$ir, m$ir)
    expect_false(.Call(zpx_model_valid, r$state$ptr))
    after <- hb$run(hb_eval_model(r), case$par, case$y)
    expect_identical(after, before)
    expect_true(.Call(zpx_model_valid, r$state$ptr))
    expect_identical(pmx_ode_tier(r), tier)
  }
})

test_that("a truncated IR is refused before anything is lowered", {
  ir <- hb_one_cpt()$ir
  # Collected first and after, so a finalizer from an earlier test running
  # in between cannot move the count.
  gc()
  live <- .Call(zpx_live_models)
  for (n in c(0L, 1L, 20L, length(ir$code) - 1L)) {
    bad <- new_ir(ir$code[seq_len(n)], ir$consts)
    expect_error(new_model(bad), class = "pmx_ode_error_ir")
  }
  gc()
  expect_identical(.Call(zpx_live_models), live)
})

test_that("an altered IR is refused, naming where", {
  ir <- hb_one_cpt()$ir
  rhs_start <- 7L + 15L + 9L * 5L
  cases <- list(
    magic = function(x) {
      x$code[1] <- 0L
      x
    },
    version = function(x) {
      x$code[2] <- 2L
      x
    },
    count = function(x) {
      x$code[3] <- -1L
      x
    },
    consts = function(x) {
      x$consts <- c(x$consts, 1)
      x
    },
    # rhs's first instruction is loadw r0 0; slot 99 is past nwork.
    operand = function(x) {
      x$code[rhs_start + 2L] <- 99L
      x
    },
    opcode = function(x) {
      x$code[rhs_start] <- 99L
      x
    }
  )
  for (k in names(cases)) {
    bad <- cases[[k]](ir)
    err <- expect_error(new_model(bad), class = "pmx_ode_error_ir")
    if (k %in% c("operand", "opcode")) {
      expect_identical(err$block, "rhs")
      expect_identical(err$insn, 1L)
    }
  }
})

test_that("a restored IR that was altered on disk is refused on use", {
  m <- new_model(hb_one_cpt()$ir)
  r <- save_restore(m)
  r$ir$code[1] <- 1L
  expect_error(model_eval(r, "pre", par = c(0, 0, 0), y = c(0, 0)),
               class = "pmx_ode_error_ir")
})
