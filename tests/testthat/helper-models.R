# Hand-built IR models (roadmap Stage 2): written in the IR assembler before
# the compiler exists, each with the same arithmetic written in R as its
# expected answer. `cases` are inputs; `expect(case)` gives what each block
# must return for them.

# A one-compartment oral model across three blocks: `pre` computes the
# individual parameters into work, `rhs` the derivatives from work, `post`
# the concentration.
hb_one_cpt <- function() {
  list(
    ir = ir_asm(list(
      pre = "
        loadp r0 0
        call1 r1 exp r0
        storew 0 r1
        loadp r0 1
        call1 r1 exp r0
        storew 1 r1
        loadp r0 2
        call1 r1 exp r0
        storew 2 r1
      ",
      rhs = "
        loadw r0 0
        loadw r1 1
        loadw r2 2
        loady r3 0
        loady r4 1
        mul r5 r0 r3
        neg r6 r5
        storeo 0 r6
        div r7 r1 r2
        mul r8 r7 r4
        sub r9 r5 r8
        storeo 1 r9
      ",
      post = "
        loady r0 1
        loadw r1 2
        div r2 r0 r1
        storeo 0 r2
      "
    ), ny = 2, npar = 3, nwork = 3, nout = c(rhs = 2, post = 1)),
    run = function(eval, par, y) {
      pre <- eval("pre", par = par, y = y)
      rhs <- eval("rhs", par = par, y = y, work = pre$work)
      post <- eval("post", par = par, y = y, work = pre$work)
      list(work = pre$work, rhs = rhs$out, post = post$out)
    },
    expect = function(par, y) {
      ka <- exp(par[1]); cl <- exp(par[2]); v <- exp(par[3])
      list(work = c(ka, cl, v),
           rhs = c(-(ka * y[1]), ka * y[1] - cl / v * y[2]),
           post = y[2] / v)
    },
    cases = list(
      list(par = c(log(0.8), log(2), log(20)), y = c(100, 5)),
      list(par = c(0, 0, 0), y = c(0, 0)),
      list(par = c(-3, 2, 5), y = c(1e-300, 1e300))
    )
  )
}

# Branches and logic: out = if (a < b && !(a == 0)) a * b else if (a >= b
# || a != a) a - b else t; then the six comparisons as 0/1. NaN follows C
# (design §6.3 and the IR's comparison rule), so `a != a` is the NaN test.
hb_branches <- function() {
  list(
    ir = ir_asm(list(rhs = "
      loadp r0 0
      loadp r1 1
      loadt r2
      const r3 0
      lt r4 r0 r1
      eq r5 r0 r3
      not r6 r5
      and r7 r4 r6
      jmpz r7 second
      mul r8 r0 r1
      storeo 0 r8
      jmp compare
      second:
      ge r9 r0 r1
      ne r10 r0 r0
      or r11 r9 r10
      jmpz r11 third
      sub r8 r0 r1
      storeo 0 r8
      jmp compare
      third:
      mov r8 r2
      storeo 0 r8
      compare:
      lt r12 r0 r1
      storeo 1 r12
      le r12 r0 r1
      storeo 2 r12
      gt r12 r0 r1
      storeo 3 r12
      ge r12 r0 r1
      storeo 4 r12
      eq r12 r0 r1
      storeo 5 r12
      ne r12 r0 r1
      storeo 6 r12
    "), npar = 2, nout = c(rhs = 7)),
    expect = function(par, t) {
      a <- par[1]; b <- par[2]
      c_lt <- !is.na(a < b) && a < b
      first <- if (c_lt && !(!is.na(a == 0) && a == 0)) a * b else
        if ((!is.na(a >= b) && a >= b) || is.nan(a)) a - b else t
      cmp <- function(x) as.numeric(!is.na(x) && x)
      c(first, cmp(a < b), cmp(a <= b), cmp(a > b), cmp(a >= b), cmp(a == b),
        as.numeric(is.na(a != b) || a != b))
    },
    cases = list(c(1, 2), c(0, 2), c(3, 2), c(2, 2), c(NaN, 1), c(1, NaN),
                 c(-Inf, Inf), c(-0, 0))
  )
}

# One instruction of each kind not covered above, into out: const, call2,
# neg on -0, div by zero, and a work round trip.
hb_ops <- function() {
  list(
    ir = ir_asm(list(rhs = "
      loadp r0 0
      loadp r1 1
      const r2 -0
      const r3 NaN
      call2 r4 pow r0 r1
      storeo 0 r4
      call2 r4 min r0 r1
      storeo 1 r4
      call2 r4 max r0 r1
      storeo 2 r4
      neg r4 r2
      storeo 3 r4
      div r4 r0 r2
      storeo 4 r4
      add r4 r0 r3
      storeo 5 r4
      storew 0 r1
      loadw r5 0
      storeo 6 r5
    "), npar = 2, nwork = 1, nout = c(rhs = 7)),
    expect = function(par) {
      # -0 made at run time: R's byte compiler folds a literal -0 to 0.
      negzero <- as.numeric("-0")
      c(par[1]^par[2], min(par), max(par), 0, par[1] / negzero, NaN, par[2])
    },
    cases = list(c(2, 10), c(-8, 1 / 3), c(0, 0), c(1, NaN), c(NaN, 0))
  )
}

hb_eval_model <- function(m) {
  function(block, ...) model_eval(m, block, ...)
}

hb_eval_ref <- function(ir) {
  function(block, ...) ir_eval(ir, block, ...)
}

# Each tier that runs here, by name.
available_tiers <- function() {
  t <- zupmxode_tiers()
  names(t)[t]
}
