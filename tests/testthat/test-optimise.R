# The optimisations (design §4.1): folding, common subexpressions, dead
# code and ivl hoisting each change how much code runs, never a value. The
# oracle is the same model compiled with the optimisation off.


test_that("each optimisation off gives the same doubles on every corpus model", {
  skip_under_gctorture()
  texts <- corpus_model_texts()
  offs <- list(list(fold = FALSE), list(cse = FALSE), list(hoist = FALSE),
               list(dce = FALSE),
               list(fold = FALSE, cse = FALSE, hoist = FALSE, dce = FALSE))
  for (name in names(texts)) {
    m <- tryCatch(pmx_ode_model(texts[[name]]), pmx_ode_error_compile = function(e) NULL)
    if (is.null(m)) next
    for (seed in 1:2) {
      inputs <- model_inputs(m, seed)
      want <- blocks_bits(m, inputs)
      for (o in offs) {
        got <- blocks_bits(variant(texts[[name]], o), inputs)
        expect_identical(got, want, label = paste(name, paste(names(o), collapse = "+")))
      }
    }
  }
})

test_that("the optimisations do what they say", {
  text <- "
    cl <- exp(tcl) * (2 * 3)
    v <- exp(tv)
    k <- cl / v
    unused <- exp(tcl) + 99
    d/dt(central) <- -k * central - cl / v * central
    cp <- central / v
  "
  code_len <- function(m, block) {
    d <- m$ir$code[7:21]
    d[(match(block, ir_blocks) - 1L) * 3L + 3L]
  }
  on <- variant(text, list())
  # Hoisting moves the parameter algebra out of rhs into ivl.
  expect_identical(on$hoisted, c("cl", "v", "k", "unused"))
  expect_gt(code_len(on, "ivl"), 0L)
  off <- variant(text, list(hoist = FALSE))
  expect_identical(code_len(off, "ivl"), 0L)
  expect_lt(code_len(on, "rhs"), code_len(off, "rhs"))
  # CSE: cl / v is computed once in rhs with hoisting off.
  no_cse <- variant(text, list(hoist = FALSE, cse = FALSE))
  expect_lt(code_len(off, "rhs"), code_len(no_cse, "rhs"))
  # Folding: (2 * 3) is one constant.
  no_fold <- variant(text, list(fold = FALSE))
  expect_lt(code_len(on, "ivl"), code_len(no_fold, "ivl"))
  # Dead code: `unused` is not computed for rhs.
  no_dce <- variant(text, list(hoist = FALSE, dce = FALSE))
  expect_lt(code_len(off, "rhs"), code_len(no_dce, "rhs"))
})

test_that("a condition that folds to a constant keeps only its branch", {
  a <- variant("if (1 > 2) k <- exp(p) else k <- 2\nd/dt(x) <- -k * x", list())
  # No call to exp() and no jump is left; the else branch's value is.
  ops <- matrix(a$ir$code[-(1:21)], nrow = 5)[1, ]
  expect_false(any(ops %in% c(20L, 40L, 41L)))
  expect_identical(run_blocks(a, 0, c(x = 3), c(p = 1))$rhs[["x"]], -6)
  # NaN is true, as in C: the then branch.
  m <- variant("if (0 / 0) k <- 1 else k <- 2\nd/dt(x) <- -k * x", list())
  expect_identical(run_blocks(m, 0, c(x = 1), double())$rhs[["x"]], -1)
})
