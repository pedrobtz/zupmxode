# The parallel region allocates nothing on the R heap (design §4.4): the
# per-subject solve is src/core/solve.c, which builds without R's headers
# (tests/c), and R's own allocations in a solve do not grow with the
# number of subjects, whatever the region does with them.

test_that("a solve's R allocations do not grow with its subjects", {
  skip_if_not(capabilities("profmem"), "R built without memory profiling")
  m <- pmx_ode_method(pmx_ode_model("d/dt(central) <- -k * wt / 70 * central"), "lsoda")
  # Both sizes past R's small-vector pools (a few subjects allocate from
  # them, a different count): 200 and 2000 subjects.
  small <- population_events(200)
  large <- population_events(2000)
  # Warm up both sizes, and keep R's byte compiler out of the count: the
  # first calls compile the solve's R functions, a few thousand
  # allocations that have nothing to do with the subjects.
  old_jit <- compiler::enableJIT(0)
  on.exit(compiler::enableJIT(old_jit), add = TRUE)
  for (i in 1:3) {
    pmx_ode_solve(m, small, c(k = 0.1))
    pmx_ode_solve(m, large, c(k = 0.1))
  }
  a <- profmem_count(pmx_ode_solve(m, small, c(k = 0.1)))
  b <- profmem_count(pmx_ode_solve(m, large, c(k = 0.1)))
  # The same calls allocate the same objects, larger for more subjects;
  # one allocation per subject would put b 1800 above a.
  expect_lt(abs(b - a), 10)
})

test_that("the per-subject solve is core code, built without R", {
  core <- test_path("..", "..", "src", "core", "solve.c")
  skip_if_not(file.exists(core), "the source tree is not here (an installed check)")
  src <- readLines(core)
  expect_false(any(grepl("#include *[<\"](R|Rinternals|R_ext/)", src)))
})
