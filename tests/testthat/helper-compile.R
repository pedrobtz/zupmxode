# Running compiled models block by block, and the R evaluation of the same
# statements that is their oracle (roadmap Stage 3).

# Every block of m at one time, state and parameter vector, as the solvers
# will run them: pre, then ivl into work, then rhs and post from that work.
run_blocks <- function(m, t, y, par) {
  pre <- model_eval(m, "pre", t = t, y = y, par = par)
  ivl <- model_eval(m, "ivl", t = t, y = y, par = par)
  rhs <- model_eval(m, "rhs", t = t, y = y, par = par, work = ivl$work)
  post <- model_eval(m, "post", t = t, y = y, par = par, work = ivl$work)
  ny <- length(m$states)
  groups <- matrix(pre$out, nrow = ny, dimnames = list(m$states, pre_groups))
  list(pre = groups, rhs = stats::setNames(rhs$out, m$states),
       post = stats::setNames(post$out, m$captures))
}

# The model's statements evaluated by R itself, in an environment holding
# the states, parameters and time: what each derivative, special and
# capture is. The math table's names that base R lacks are R's own
# equivalents. linCmt() has no R meaning, so such models are not run here.
r_eval_model <- function(text, t, y, par, states, captures) {
  env <- list2env(c(as.list(y), as.list(par), list(t = t, time = t)),
                  parent = baseenv())
  env$expit <- function(x, low = 0, high = 1) low + (high - low) * stats::plogis(x)
  env$logit <- function(x, low = 0, high = 1) stats::qlogis((x - low) / (high - low))
  env$probit <- stats::qnorm
  env$probitInv <- stats::pnorm
  out <- list(rhs = stats::setNames(numeric(length(states)), states),
              pre = matrix(rep(c(0, 1, 0, 0, 0), each = length(states)),
                           nrow = length(states),
                           dimnames = list(states, pre_groups)))
  key <- c(alag = "lag", lag = "lag", f = "f", F = "f", rate = "rate", dur = "dur")
  run <- function(e) {
    if (is_call_to(e, "cmt")) {
      return(invisible())
    } else if (is_call_to(e, "{")) {
      for (x in as.list(e)[-1]) run(x)
    } else if (is_call_to(e, "if")) {
      if (eval(e[[2]], env)) run(e[[3]]) else if (length(e) == 4) run(e[[4]])
    } else {
      lhs <- e[[2]]
      val <- eval(e[[3]], env)
      if (is.name(lhs)) {
        assign(as.character(lhs), val, envir = env)
      } else if (is_call_to(lhs, "/")) {
        out$rhs[[as.character(lhs[[3]][[2]])]] <<- val
      } else if (identical(lhs[[2]], 0) || identical(lhs[[2]], 0L)) {
        out$pre[as.character(lhs[[1]]), "init"] <<- val
      } else {
        out$pre[as.character(lhs[[2]]), key[[as.character(lhs[[1]])]]] <<- val
      }
    }
  }
  for (e in as.list(parse(text = text))) run(e)
  out$post <- vapply(captures, function(n) get(n, envir = env), 0)
  out
}

# Inputs for a compiled model: states and parameters drawn so that every
# value is positive and moderate, as PK parameters are.
model_inputs <- function(m, seed) {
  withr::with_seed(seed, list(
    t = stats::runif(1, 0, 24),
    y = stats::setNames(stats::runif(length(m$states), 0.1, 50), m$states),
    par = stats::setNames(stats::runif(length(m$params), 0.1, 1.5), m$params)
  ))
}

# The oracle corpus's model texts, by name.
corpus_model_texts <- function() {
  dir <- test_path("fixtures", "oracle", "models")
  fs <- list.files(dir, pattern = "[.]R$")
  stats::setNames(lapply(file.path(dir, fs), function(f) paste(readLines(f), collapse = "\n")),
                  sub("[.]R$", "", fs))
}

uses_lincmt <- function(text) grepl("linCmt", text, fixed = TRUE)

# ---- Shared by the Stage 3 test files (the suite runs shuffled, so nothing a
# test uses may live at a test file's top level) ----

# From test-classify.R.
# The classifier (design §4.1): closed-form for linCmt() with no ODE;
# linear when every derivative is affine in the states with coefficients
# constant within an interval; general otherwise. The labels below were
# written by reading each model, not by running the classifier.

hand_labels <- c(
  emax_effect = "linear", irm1 = "general", pk1_iv = "linear",
  pk1_iv_lag = "linear", pk1_iv_ratedur = "linear", pk1_iv_wt = "linear",
  pk1_mm = "general", pk1_oral = "linear", pk1_oral_lagf = "linear",
  pk1_oral_lin = "closed-form", pk2_iv_lin = "closed-form",
  pk2_oral = "linear", pk3_iv_lin = "closed-form",
  lib_PK_1cmt = "closed-form", lib_PK_1cmt_des = "linear",
  lib_PK_2cmt = "closed-form", lib_PK_2cmt_des = "linear",
  lib_PK_2cmt_mAb_Davda_2014 = "closed-form", lib_PK_2cmt_no_depot = "linear",
  # clearance changes with time
  lib_PK_2cmt_tdcl_des = "general",
  lib_PK_3cmt = "closed-form", lib_PK_3cmt_des = "linear",
  lib_PK_double_sim_01 = "linear", lib_PK_double_sim_10 = "linear",
  lib_PK_double_sim_11 = "linear",
  # brain compartments driven linearly by plasma
  lib_Grimm_2023_gantenerumab = "linear",
  # Michaelis-Menten elimination
  lib_Grimm_2023_trontinemab = "general",
  lib_Kovalenko_2020_dupilumab = "general",
  lib_indirect_0cpt_transitEx = "general",
  lib_indirect_1cpt_inhi_kin = "general", lib_indirect_1cpt_inhi_kin_CLV = "general",
  lib_indirect_1cpt_inhi_kin_r0rmaxcrmax = "general",
  lib_indirect_1cpt_inhi_kout_CLV = "general",
  lib_indirect_1cpt_inhi_kout_r0rmaxcrmax = "general",
  lib_indirect_1cpt_stim_kout_CLV = "general",
  lib_indirect_prec_1cpt_inhi_CLV = "general",
  lib_indirect_prec_1cpt_stim_CLV = "general",
  lib_tgi_no_sat_Koch = "general",
  # constant growth: affine
  lib_tgi_no_sat_linear = "linear",
  lib_tgi_no_sat_powerLaw = "general", lib_tgi_sat_Gompertz = "general",
  lib_tgi_sat_VonBertalanffy = "general", lib_tgi_sat_genVonBertalanffy = "general"
)

# From test-compile-errors.R.
# Compile errors (design §8): a construct outside the subset is an error of
# class pmx_ode_error_compile naming the construct, its position and the
# nearest supported alternative.

expect_compile_error <- function(text, construct, alternative = NULL) {
  err <- expect_error(pmx_ode_model(text), class = "pmx_ode_error_compile")
  expect_identical(err$construct, construct, label = text)
  expect_false(is.na(err$line), label = text)
  expect_false(is.na(err$column), label = text)
  if (!is.null(alternative)) expect_match(err$alternative, alternative, fixed = TRUE)
  expect_match(conditionMessage(err), paste0("line ", err$line), fixed = TRUE)
  invisible(err)
}

# From test-language.R.
# The language (design §6.1): every construct compiles, and what it computes
# is what R computes for the same statements on the same values. R is the
# oracle because the language is R syntax; rxode2 is the second one for
# the promise that what compiles here compiles there.

construct_cases <- function() {
  c("a + b", "a - b", "a * b", "a / b", "a ^ b", "a ** b", "-a", "+a",
    "(a + b) * a", "exp(a)", "log(a)", "sqrt(a)", "abs(-a)", "floor(a * 10)",
    "ceiling(a * 10)", "sin(a)", "cos(a)", "tan(a)", "gamma(a)", "lgamma(a)",
    "logit(a / 2)", "expit(a)", "probit(a / 2)", "probitInv(a)", "min(a, b)",
    "max(a, b)", "min(a, b, 1)", "max(1, a, b)", "a < b", "a <= b", "a > b",
    "a >= b", "a == b", "a != b", "a < b && b > 1", "a > b || b > 1",
    "a < b & b > 1", "a > b | b > 1", "!(a < b)", "TRUE", "FALSE", "pi",
    "t * a", "time * a", "logit(a, 0, 2)", "expit(a, 1, 3)", "x * a",
    "exp(a) * exp(a) + exp(a)", "1 / 0", "-1 / 0", "2 ^ 0.5 * a")
}

# From test-language.R.

expect_matches_r <- function(text, inputs, label, tolerance = 1e-14) {
  m <- pmx_ode_model(text)
  got <- run_blocks(m, inputs$t, inputs$y[m$states], inputs$par[m$params])
  want <- r_eval_model(text, inputs$t, inputs$y[m$states],
                       inputs$par[m$params], m$states, m$captures)
  for (part in c("rhs", "pre", "post")) {
    expect_equal(unname(got[[part]]), unname(want[[part]]),
                 tolerance = tolerance, label = paste(label, part))
  }
}

# From test-optimise.R.
# The optimisations (design §4.1): folding, common subexpressions, dead
# code and ivl hoisting each change how much code runs, never a value. The
# oracle is the same model compiled with the optimisation off.

variant <- function(text, opts) {
  new_compiled_model(compile_model(text, opts))
}

# From test-optimise.R.

blocks_bits <- function(m, inputs) {
  b <- run_blocks(m, inputs$t, inputs$y[m$states], inputs$par[m$params])
  lapply(b, function(x) writeBin(as.vector(x), raw()))
}

# The sweeps over the whole corpus run R code: the compiler, and R's own
# evaluation as their oracle. gctorture is for the C entries, which
# test-tiers.R, test-model-object.R and test-ir-serialise.R exercise under
# it; the sweeps run on every other leg, and at step 500 they alone would
# take past gctorture's 10-minute budget.
skip_under_gctorture <- function() {
  testthat::skip_if(identical(Sys.getenv("R_ACTIONS_WORKFLOW"), "gctorture"),
                    "a corpus sweep of R code, run on every other leg")
}
