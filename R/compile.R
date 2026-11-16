# The compiler (design section 4.1, section 6.1): a model in the rxode2 subset, read by R's
# parser, becomes IR in five blocks. Statements are checked against the
# subset, names resolved into states, parameters, locals and captures, the
# system classified, and each block generated from only the statements it
# needs, with constants folded, common subexpressions merged and
# state-independent statements of `rhs` hoisted into `ivl`.

# ---- Reading the model ------------------------------------------------------

# The model's top-level statements, each with its line and column when the
# source is known.
model_statements <- function(model, env) {
  if (is.character(model)) {
    exprs <- parse(text = model, keep.source = TRUE)
    refs <- attr(exprs, "srcref")
    exprs <- as.list(exprs)
  } else {
    if (is.function(model)) model <- body(model)
    if (is.call(model) && identical(model[[1]], as.name("{"))) {
      refs <- attr(model, "srcref")
      exprs <- as.list(model)[-1]
      refs <- if (is.null(refs)) NULL else refs[-1]
    } else if (is.expression(model)) {
      refs <- attr(model, "srcref")
      exprs <- as.list(model)
    } else {
      exprs <- list(model)
      refs <- NULL
    }
  }
  lapply(seq_along(exprs), function(i) {
    ref <- if (length(refs) >= i) refs[[i]] else NULL
    list(expr = exprs[[i]],
         line = if (is.null(ref)) NA_integer_ else as.integer(ref[1]),
         col = if (is.null(ref)) NA_integer_ else as.integer(ref[2]),
         index = i)
  })
}

# ---- Compile errors ---------------------------------------------------------

# rxode2 constructs outside the subset (design section 6.1), with the rxode2
# feature each belongs to and the nearest supported alternative.
rxode2_only <- list(
  ifelse = c("rxode2's ifelse()", "if (...) ... else ..."),
  `~` = c("a residual error line", "residual error given to pmx_ode_sim(sigma = )"),
  `for` = c("a loop", "the steps written out"),
  `while` = c("a loop", "the steps written out"),
  `repeat` = c("a loop", "the steps written out"),
  transit = c("rxode2's transit compartments", "the compartments as d/dt() equations"),
  mtime = c("a model event time", NA),
  sum = c("rxode2's sum()", "the terms written out with +"),
  prod = c("rxode2's prod()", "the factors written out with *"),
  log10 = c("a function outside the math table", "log(x) / log(10)"),
  log2 = c("a function outside the math table", "log(x) / log(2)"),
  log1p = c("a function outside the math table", "log(1 + x)"),
  expm1 = c("a function outside the math table", "exp(x) - 1"),
  rxnorm = c("a random draw inside the model", "draws made by pmx_ode_sim() (design D7)"),
  rnorm = c("a random draw inside the model", "draws made by pmx_ode_sim() (design D7)"),
  rxbinom = c("a random draw inside the model", "draws made by pmx_ode_sim() (design D7)"),
  runif = c("a random draw inside the model", "draws made by pmx_ode_sim() (design D7)"),
  `<<-` = c("a global assignment", "<-"),
  `->` = c("a right assignment", "<-"),
  `[` = c("indexing", NA),
  `[[` = c("indexing", NA),
  `$` = c("indexing", NA)
)

abort_compile <- function(construct, st, what, alternative = NA, call = NULL) {
  where <- if (is.na(st$line)) {
    paste0("statement ", st$index)
  } else {
    paste0("line ", st$line, ", column ", st$col)
  }
  pmx_abort(
    "pmx_ode_error_compile",
    paste0("Unsupported construct `", construct, "` at ", where, ": ", what,
           if (!is.na(alternative)) paste0(". Use ", alternative, " instead"),
           "."),
    construct = construct, line = st$line, column = st$col,
    statement = st$index, alternative = alternative, call = call
  )
}

unsupported <- function(construct, st) {
  hit <- rxode2_only[[construct]]
  if (!is.null(hit)) {
    abort_compile(construct, st,
                  paste0(hit[1], ", an rxode2 feature outside zupmxode's subset"),
                  hit[2])
  }
  near <- agrep(construct, c(names(ir_math), "linCmt"), max.distance = 0.2,
                value = TRUE)
  abort_compile(construct, st, "not in zupmxode's model language (design section 6.1)",
                if (length(near)) paste0(near[1], "()") else NA)
}

# ---- Statements ------------------------------------------------------------

special_fns <- c(alag = "lag", lag = "lag", f = "f", F = "f", rate = "rate",
                 dur = "dur")

is_call_to <- function(e, name) is.call(e) && identical(e[[1]], as.name(name))

# One parsed statement as a node: assign, deriv, special, cmt or if.
parse_stmt <- function(e, st) {
  if (is_call_to(e, "{")) {
    return(lapply(as.list(e)[-1], parse_stmt, st = st))
  }
  if (is_call_to(e, "if")) {
    then <- parse_block(e[[3]], st)
    els <- if (length(e) == 4) parse_block(e[[4]], st) else list()
    return(list(list(kind = "if", cond = check_expr(e[[2]], st), then = then,
                     els = els, st = st)))
  }
  if (is_call_to(e, "<-") || is_call_to(e, "=")) {
    lhs <- e[[2]]
    rhs <- e[[3]]
    if (is.name(lhs)) {
      name <- as.character(lhs)
      if (name %in% c("t", "time", "pi")) {
        abort_compile(paste0(name, " <-"), st, "an assignment to a reserved name",
                      "another name")
      }
      if (is_call_to(rhs, "linCmt")) {
        if (length(rhs) > 1) {
          abort_compile("linCmt(...)", st, "linCmt() with arguments", "linCmt()")
        }
        return(list(list(kind = "lincmt", name = name, st = st)))
      }
      return(list(list(kind = "assign", name = name, expr = check_expr(rhs, st),
                       st = st)))
    }
    # d/dt(x) <- e
    if (is_call_to(lhs, "/") && identical(lhs[[2]], as.name("d")) &&
          is_call_to(lhs[[3]], "dt") && length(lhs[[3]]) == 2 &&
          is.name(lhs[[3]][[2]])) {
      return(list(list(kind = "deriv", state = as.character(lhs[[3]][[2]]),
                       expr = check_expr(rhs, st), st = st)))
    }
    if (is.call(lhs) && is.name(lhs[[1]]) && length(lhs) == 2) {
      fn <- as.character(lhs[[1]])
      arg <- lhs[[2]]
      # x(0) <- e: an initial condition.
      if (identical(arg, 0) || identical(arg, 0L)) {
        return(list(list(kind = "special", what = "init", state = fn,
                         expr = check_expr(rhs, st), st = st)))
      }
      if (fn %in% names(special_fns) && is.name(arg)) {
        return(list(list(kind = "special", what = special_fns[[fn]],
                         state = as.character(arg), expr = check_expr(rhs, st),
                         st = st)))
      }
    }
    abort_compile(paste(deparse(lhs), "<-"), st, "an assignment target outside the subset",
                  "a name, d/dt(x), x(0), alag(x), f(x), rate(x) or dur(x)")
  }
  if (is_call_to(e, "cmt") && length(e) == 2 && is.name(e[[2]])) {
    return(list(list(kind = "cmt", state = as.character(e[[2]]), st = st)))
  }
  if (is.call(e) && is.name(e[[1]])) unsupported(as.character(e[[1]]), st)
  abort_compile(paste(deparse(e), collapse = " "), st,
                "a statement that is not an assignment",
                "an assignment with <- or =")
}

parse_block <- function(e, st) {
  out <- parse_stmt(e, st)
  if (!is.null(out$kind)) list(out) else unlist_nodes(out)
}

unlist_nodes <- function(x) {
  out <- list()
  for (n in x) {
    if (!is.null(n$kind)) out[[length(out) + 1L]] <- n else out <- c(out, unlist_nodes(n))
  }
  out
}

# ---- Expressions -----------------------------------------------------------

binary_ops <- c(`+` = "add", `-` = "sub", `*` = "mul", `/` = "div",
                `^` = "pow", `<` = "lt", `<=` = "le", `>` = "gt", `>=` = "ge",
                `==` = "eq", `!=` = "ne", `&&` = "and", `&` = "and",
                `||` = "or", `|` = "or")

# An expression as a tree of list(op, args) nodes, checked against the
# subset: numbers, names, the operators above, unary - and !, and the math
# table's functions.
check_expr <- function(e, st) {
  if (is.numeric(e) && length(e) == 1L) return(list(op = "num", value = as.double(e)))
  if (is.logical(e) && length(e) == 1L && !is.na(e)) {
    return(list(op = "num", value = as.double(e)))
  }
  if (is.name(e)) {
    n <- as.character(e)
    if (n == "pi") return(list(op = "num", value = pi))
    if (n == "Inf") return(list(op = "num", value = Inf))
    if (n == "NaN") return(list(op = "num", value = NaN))
    if (n %in% c("t", "time")) return(list(op = "time"))
    return(list(op = "name", name = n))
  }
  if (!is.call(e)) abort_compile(paste(deparse(e), collapse = " "), st, "a value that is not a number or a name")
  if (!is.name(e[[1]])) abort_compile(paste(deparse(e[[1]]), collapse = " "), st, "a call of something that is not a function name")
  fn <- as.character(e[[1]])
  args <- as.list(e)[-1]
  if (fn == "(" && length(args) == 1L) return(check_expr(args[[1]], st))
  if (fn %in% c("-", "+") && length(args) == 1L) {
    a <- check_expr(args[[1]], st)
    return(if (fn == "-") list(op = "neg", args = list(a)) else a)
  }
  if (fn == "!" && length(args) == 1L) {
    return(list(op = "not", args = list(check_expr(args[[1]], st))))
  }
  if (fn %in% names(binary_ops) && length(args) == 2L) {
    return(list(op = binary_ops[[fn]], args = lapply(args, check_expr, st = st)))
  }
  if (fn %in% c("min", "max") && length(args) >= 2L) {
    a <- lapply(args, check_expr, st = st)
    out <- a[[1]]
    for (b in a[-1]) out <- list(op = "call", fn = fn, args = list(out, b))
    return(out)
  }
  # rxode2's bounded logit(x, low, high) and expit(x, low, high), written
  # through the one-argument functions: logit((x - low) / (high - low)) and
  # low + (high - low) * expit(x).
  if (fn %in% c("logit", "expit") && length(args) %in% c(2L, 3L)) {
    a <- lapply(args, check_expr, st = st)
    low <- a[[2]]
    high <- if (length(a) == 3L) a[[3]] else list(op = "num", value = 1)
    width <- list(op = "sub", args = list(high, low))
    if (fn == "logit") {
      inner <- list(op = "div", args = list(list(op = "sub", args = list(a[[1]], low)), width))
      return(list(op = "call", fn = "logit", args = list(inner)))
    }
    scaled <- list(op = "mul", args = list(width, list(op = "call", fn = "expit", args = list(a[[1]]))))
    return(list(op = "add", args = list(low, scaled)))
  }
  if (fn %in% names(ir_math) && ir_math[[fn]] < 32L) {
    if (length(args) != 1L) {
      abort_compile(paste0(fn, "()"), st, paste0(fn, "() with ", length(args), " arguments"),
                    paste0(fn, "(x)"))
    }
    return(list(op = "call", fn = fn, args = list(check_expr(args[[1]], st))))
  }
  if (fn == "linCmt") {
    abort_compile("linCmt()", st, "linCmt() inside an expression", "name <- linCmt()")
  }
  unsupported(fn, st)
}

expr_names <- function(x) {
  if (x$op == "name") return(x$name)
  if (is.null(x$args)) return(character())
  unique(unlist(lapply(x$args, expr_names)))
}

uses_time <- function(x) {
  if (x$op == "time") return(TRUE)
  if (is.null(x$args)) return(FALSE)
  any(vapply(x$args, uses_time, NA))
}

# ---- linCmt() --------------------------------------------------------------

# The parameterisations linCmt() recognises (design section 4.3), by the names the
# model assigns, as rxode2 spells them (case is ignored, as there).
lincmt_names <- list(
  cl = c("cl"), v = c("v", "vc", "v1"), ka = c("ka"), q = c("q", "q1", "cld"),
  vp = c("vp", "v2", "vt", "vp1"), q2 = c("q2"), vp2 = c("vp2", "v3"),
  k = c("k", "kel", "ke"), k12 = c("k12"), k21 = c("k21"), k13 = c("k13"),
  k31 = c("k31")
)

lincmt_resolve <- function(assigned, st) {
  low <- tolower(assigned)
  pick <- function(key) {
    hit <- assigned[low %in% lincmt_names[[key]]]
    if (length(hit) > 1L) {
      abort_compile("linCmt()", st, paste0("more than one name for ", key, ": ",
                                           paste(hit, collapse = ", ")))
    }
    if (length(hit)) hit else NA_character_
  }
  p <- vapply(names(lincmt_names), pick, "")
  if (is.na(p[["v"]])) {
    abort_compile("linCmt()", st, "linCmt() with no central volume assigned before it",
                  "v (or vc, v1) assigned before linCmt()")
  }
  micro <- !is.na(p[["k"]])
  if (!micro && is.na(p[["cl"]])) {
    abort_compile("linCmt()", st, "linCmt() with neither a clearance nor an elimination rate",
                  "cl or k assigned before linCmt()")
  }
  # (!is.na(.)) in parentheses: `!` binds more loosely than `+`.
  ncmt <- if (micro) {
    1L + (!is.na(p[["k12"]])) + (!is.na(p[["k13"]]))
  } else {
    1L + (!is.na(p[["q"]])) + (!is.na(p[["q2"]]))
  }
  list(params = p, micro = micro, ncmt = ncmt, depot = !is.na(p[["ka"]]))
}

# linCmt() as the ODE system it solves, as statements: rxode2's compartment
# names, so a model's doses and observations name the same compartments the
# closed form will (Stage 5) and the general solvers do now.
lincmt_statements <- function(lc, name, st) {
  p <- lc$params
  nm <- function(x) list(op = "name", name = x)
  op <- function(o, a, b) list(op = o, args = list(a, b))
  div <- function(a, b) op("div", nm(a), nm(b))
  k <- if (lc$micro) nm(p[["k"]]) else div(p[["cl"]], p[["v"]])
  k12 <- if (lc$micro) nm(p[["k12"]]) else if (lc$ncmt >= 2) div(p[["q"]], p[["v"]])
  k21 <- if (lc$micro) nm(p[["k21"]]) else if (lc$ncmt >= 2) div(p[["q"]], p[["vp"]])
  k13 <- if (lc$micro) nm(p[["k13"]]) else if (lc$ncmt >= 3) div(p[["q2"]], p[["v"]])
  k31 <- if (lc$micro) nm(p[["k31"]]) else if (lc$ncmt >= 3) div(p[["q2"]], p[["vp2"]])
  central <- op("mul", list(op = "neg", args = list(k)), nm("central"))
  if (lc$depot) central <- op("add", op("mul", nm(p[["ka"]]), nm("depot")), central)
  out <- list()
  if (lc$depot) {
    out[[length(out) + 1L]] <- list(kind = "deriv", state = "depot", st = st,
                                     expr = list(op = "neg", args = list(op("mul", nm(p[["ka"]]), nm("depot")))))
  }
  if (lc$ncmt >= 2) {
    central <- op("add", op("sub", central, op("mul", k12, nm("central"))),
                  op("mul", k21, nm("peripheral1")))
  }
  if (lc$ncmt >= 3) {
    central <- op("add", op("sub", central, op("mul", k13, nm("central"))),
                  op("mul", k31, nm("peripheral2")))
  }
  out[[length(out) + 1L]] <- list(kind = "deriv", state = "central", expr = central, st = st)
  if (lc$ncmt >= 2) {
    out[[length(out) + 1L]] <- list(kind = "deriv", state = "peripheral1", st = st,
                                     expr = op("sub", op("mul", k12, nm("central")),
                                               op("mul", k21, nm("peripheral1"))))
  }
  if (lc$ncmt >= 3) {
    out[[length(out) + 1L]] <- list(kind = "deriv", state = "peripheral2", st = st,
                                     expr = op("sub", op("mul", k13, nm("central")),
                                               op("mul", k31, nm("peripheral2"))))
  }
  out[[length(out) + 1L]] <- list(kind = "assign", name = name, st = st,
                                   expr = op("div", nm("central"), nm(p[["v"]])))
  out
}

# ---- Names -----------------------------------------------------------------

assigned_names <- function(nodes) {
  out <- character()
  for (n in nodes) {
    if (n$kind %in% c("assign", "lincmt")) out <- c(out, n$name)
    if (n$kind == "if") out <- c(out, assigned_names(n$then), assigned_names(n$els))
  }
  unique(out)
}

deriv_states <- function(nodes) {
  out <- character()
  for (n in nodes) {
    if (n$kind %in% c("deriv", "cmt")) out <- c(out, n$state)
    if (n$kind == "if") out <- c(out, deriv_states(n$then), deriv_states(n$els))
  }
  unique(out)
}

# Checks that every name read is a state, a parameter or a local assigned
# on every path before it; returns the parameters in order of first use.
resolve_names <- function(nodes, states, locals) {
  params <- character()
  walk <- function(nodes, defined) {
    for (n in nodes) {
      reads <- if (!is.null(n$expr)) expr_names(n$expr) else if (n$kind == "if") expr_names(n$cond) else character()
      for (r in reads) {
        if (r %in% states || r %in% defined) next
        if (r %in% locals) {
          abort_compile(r, n$st, paste0("`", r, "` may be read before it is assigned"),
                        paste0("an assignment to `", r, "` on every path before"))
        }
        if (!r %in% params) params <<- c(params, r)
      }
      if (n$kind %in% c("assign", "lincmt")) defined <- union(defined, n$name)
      if (n$kind == "if") {
        a <- walk(n$then, defined)
        b <- walk(n$els, defined)
        defined <- union(defined, intersect(a, b))
      }
      if (n$kind == "special" && !n$state %in% states) {
        abort_compile(paste0(n$what, "(", n$state, ")"), n$st,
                      paste0("`", n$state, "` is not a compartment of the model"),
                      paste0("d/dt(", n$state, ") or cmt(", n$state, ")"))
      }
    }
    defined
  }
  walk(nodes, character())
  params
}

# ---- Classification --------------------------------------------------------

# The kind of an expression in the states: 0 constant within an interval (no
# state, no time), 1 affine in the states with such coefficients, 2 anything
# else.
expr_kind <- function(x, kinds) {
  switch(x$op,
    num = 0L,
    time = 2L,
    name = if (x$name %in% names(kinds)) kinds[[x$name]] else 0L,
    add = , sub = max(vapply(x$args, expr_kind, 0L, kinds = kinds)),
    neg = expr_kind(x$args[[1]], kinds),
    mul = {
      k <- vapply(x$args, expr_kind, 0L, kinds = kinds)
      if (min(k) == 0L) max(k) else 2L
    },
    div = {
      k <- vapply(x$args, expr_kind, 0L, kinds = kinds)
      if (k[2] == 0L) k[1] else 2L
    },
    if (all(vapply(x$args, expr_kind, 0L, kinds = kinds) == 0L)) 0L else 2L
  )
}

classify <- function(nodes, states, has_lincmt) {
  kinds <- stats::setNames(rep(1L, length(states)), states)
  derivs <- 0L
  walk <- function(nodes, cond_kind) {
    for (n in nodes) {
      if (n$kind == "assign") kinds[[n$name]] <<- max(expr_kind(n$expr, kinds), cond_kind)
      if (n$kind == "deriv") {
        k <- max(expr_kind(n$expr, kinds), cond_kind)
        derivs <<- max(derivs, k)
      }
      if (n$kind == "if") {
        ck <- max(cond_kind, if (expr_kind(n$cond, kinds) == 0L) 0L else 2L)
        walk(n$then, ck)
        walk(n$els, ck)
      }
    }
  }
  walk(nodes, 0L)
  user_derivs <- any(vapply(nodes, function(n) n$kind == "deriv" && !isTRUE(n$st$generated), NA))
  if (has_lincmt && !user_derivs) return("closed-form")
  if (derivs <= 1L) "linear" else "general"
}

# ---- Blocks ----------------------------------------------------------------

# Every derivative and special as an assignment to a variable no model can
# name (".d.x", ".lag.x", ...), so one code generator serves all of them.
as_assignments <- function(nodes) {
  lapply(nodes, function(n) {
    if (n$kind == "deriv") return(list(kind = "assign", name = paste0(".d.", n$state), expr = n$expr, st = n$st))
    if (n$kind == "special") return(list(kind = "assign", name = paste0(".", n$what, ".", n$state), expr = n$expr, st = n$st))
    if (n$kind == "if") {
      n$then <- as_assignments(n$then)
      n$els <- as_assignments(n$els)
    }
    n
  })
}

node_reads <- function(n) {
  if (n$kind == "assign") return(expr_names(n$expr))
  if (n$kind == "if") return(unique(c(expr_names(n$cond), unlist(lapply(c(n$then, n$els), node_reads)))))
  character()
}

node_writes <- function(n) {
  if (n$kind == "assign") return(n$name)
  if (n$kind == "if") return(unique(unlist(lapply(c(n$then, n$els), node_writes))))
  character()
}

node_time <- function(n) {
  if (n$kind == "assign") return(uses_time(n$expr))
  if (n$kind == "if") return(uses_time(n$cond) || any(vapply(c(n$then, n$els), node_time, NA)))
  FALSE
}

# The statements of `nodes` needed to compute `want`, in order (dead code
# eliminated), or every statement when dce is off.
needed <- function(nodes, want, dce = TRUE) {
  if (!dce) return(nodes)
  keep <- logical(length(nodes))
  need <- want
  for (i in rev(seq_along(nodes))) {
    w <- node_writes(nodes[[i]])
    if (any(w %in% need)) {
      keep[i] <- TRUE
      if (nodes[[i]]$kind == "assign") need <- setdiff(need, w)
      need <- union(need, node_reads(nodes[[i]]))
    }
  }
  out <- nodes[keep]
  lapply(out, function(n) {
    if (n$kind == "if") {
      n$then <- needed(n$then, union(want, unlist(lapply(out, node_reads))), dce)
      n$els <- needed(n$els, union(want, unlist(lapply(out, node_reads))), dce)
    }
    n
  })
}

# Top-level assignments, each the only assignment to its name, that read no
# state and no time, directly or through a local: what `ivl` computes once
# per interval for `rhs` and `post` (design section 4.1).
hoistable <- function(nodes, states) {
  writes <- unlist(lapply(nodes, node_writes))
  once <- names(which(table(writes) == 1L))
  varying <- states
  out <- character()
  for (n in nodes) {
    if (n$kind == "assign" && n$name %in% once && !startsWith(n$name, ".") &&
          !node_time(n) && !any(expr_names(n$expr) %in% varying)) {
      out <- c(out, n$name)
    } else {
      varying <- union(varying, node_writes(n))
    }
  }
  out
}

# ---- Code generation -------------------------------------------------------

# Generates one block as IR assembler text: `nodes` run in order; `inputs`
# maps names to where they are read from (list(kind = "y"|"p"|"w", index));
# `outputs` the variables stored to out[0..], or to work when store_work.
gen_block <- function(nodes, inputs, outputs, store_work = FALSE,
                      fold = TRUE, cse = TRUE) {
  lines <- character()
  nreg <- 0L
  emit <- function(...) lines[[length(lines) + 1L]] <<- paste(...)
  newreg <- function() {
    nreg <<- nreg + 1L
    paste0("r", nreg - 1L)
  }
  label_n <- 0L
  newlabel <- function() {
    label_n <<- label_n + 1L
    paste0("L", label_n)
  }
  memo <- new.env(parent = emptyenv())
  vals <- list()   # name -> list(reg =) or list(const =)
  # A value as a register.
  as_reg <- function(v) {
    if (!is.null(v$reg)) return(v$reg)
    key <- paste0("const:", sprintf("%a", v$const))
    if (cse && !is.null(memo[[key]])) return(memo[[key]])
    r <- newreg()
    # The bits, not a decimal: see ir_number().
    emit("const", r, ir_bits(v$const))
    if (cse) memo[[key]] <- r
    r
  }
  cached <- function(key, make) {
    if (cse && !is.null(memo[[key]])) return(list(reg = memo[[key]]))
    r <- make()
    if (cse) memo[[key]] <- r
    list(reg = r)
  }
  read_input <- function(name) {
    src <- inputs[[name]]
    op <- switch(src$kind, y = "loady", p = "loadp", w = "loadw")
    cached(paste0(op, ":", src$index), function() {
      r <- newreg()
      emit(op, r, src$index)
      r
    })
  }
  fold_value <- function(op, a, b = NULL, fn = NULL) {
    x <- a$const
    y <- if (is.null(b)) NULL else b$const
    c_cmp <- function(z) if (is.na(z)) 0 else as.double(z)
    switch(op,
      add = x + y, sub = x - y, mul = x * y, div = x / y, neg = -x,
      lt = c_cmp(x < y), le = c_cmp(x <= y), gt = c_cmp(x > y), ge = c_cmp(x >= y),
      eq = c_cmp(x == y), ne = if (is.na(x != y)) 1 else as.double(x != y),
      # C's truth: any value but 0 is true, NaN included.
      and = as.double((is.nan(x) || x != 0) && (is.nan(y) || y != 0)),
      or = as.double((is.nan(x) || x != 0) || (is.nan(y) || y != 0)),
      not = as.double(!is.nan(x) && x == 0),
      call = .Call(zpx_math_r, ir_math[[fn]], x, if (is.null(y)) double() else y),
      pow = .Call(zpx_math_r, ir_math[["pow"]], x, y)
    )
  }
  gen <- function(x) {
    switch(x$op,
      num = list(const = x$value),
      time = cached("loadt", function() {
        r <- newreg()
        emit("loadt", r)
        r
      }),
      name = {
        v <- vals[[x$name]]
        if (!is.null(v)) v else read_input(x$name)
      },
      {
        args <- lapply(x$args, gen)
        if (fold && all(vapply(args, function(a) !is.null(a$const), NA))) {
          return(list(const = fold_value(x$op, args[[1]], if (length(args) > 1) args[[2]], x$fn)))
        }
        regs <- vapply(args, as_reg, "")
        if (x$op == "call" || x$op == "pow") {
          fn <- if (x$op == "pow") "pow" else x$fn
          two <- ir_math[[fn]] >= 32L
          return(cached(paste(fn, paste(regs, collapse = ",")), function() {
            r <- newreg()
            emit(if (two) "call2" else "call1", r, fn, paste(regs, collapse = " "))
            r
          }))
        }
        cached(paste(x$op, paste(regs, collapse = ",")), function() {
          r <- newreg()
          emit(x$op, r, paste(regs, collapse = " "))
          r
        })
      }
    )
  }
  run <- function(nodes) {
    for (n in nodes) {
      if (n$kind == "assign") {
        vals[[n$name]] <<- gen(n$expr)
      } else if (n$kind == "if") {
        cond <- gen(n$cond)
        if (fold && !is.null(cond$const)) {
          # JMPZ jumps only on 0: NaN takes the then branch, as in C.
          run(if (is.nan(cond$const) || cond$const != 0) n$then else n$els)
          next
        }
        creg <- as_reg(cond)
        writes <- unique(c(unlist(lapply(n$then, node_writes)), unlist(lapply(n$els, node_writes))))
        before <- vals
        memo_before <- as.list(memo)
        merge <- stats::setNames(lapply(writes, function(w) newreg()), writes)
        branch <- function(body) {
          rm(list = ls(memo), envir = memo)
          for (k in names(memo_before)) memo[[k]] <- memo_before[[k]]
          vals <<- before
          run(body)
          for (w in writes) {
            v <- vals[[w]]
            src <- if (is.null(v)) {
              # Not assigned on this path: its value before the if, or 0 for
              # the outputs that start at 0.
              b <- before[[w]]
              if (is.null(b)) list(const = 0) else b
            } else v
            emit("mov", merge[[w]], as_reg(src))
          }
        }
        l_else <- newlabel()
        l_end <- newlabel()
        emit("jmpz", creg, l_else)
        branch(n$then)
        emit("jmp", l_end)
        emit(paste0(l_else, ":"))
        branch(n$els)
        emit(paste0(l_end, ":"))
        rm(list = ls(memo), envir = memo)
        for (k in names(memo_before)) memo[[k]] <- memo_before[[k]]
        vals <<- before
        for (w in writes) vals[[w]] <<- list(reg = merge[[w]])
      }
    }
  }
  run(nodes)
  for (i in seq_along(outputs)) {
    v <- vals[[outputs[[i]]]]
    if (is.null(v)) v <- if (outputs[[i]] %in% names(inputs)) read_input(outputs[[i]]) else list(const = 0)
    emit(if (store_work) "storew" else "storeo", i - 1L, as_reg(v))
  }
  paste(lines, collapse = "\n")
}

# ---- The compiler ----------------------------------------------------------

# pre's outputs: per state, its lag, bioavailability, rate, duration and
# initial amount, in that order of groups.
pre_groups <- c("lag", "f", "rate", "dur", "init")
pre_defaults <- c(lag = 0, f = 1, rate = 0, dur = 0, init = 0)

compile_model <- function(model, opts = list()) {
  fold <- !isFALSE(opts$fold)
  cse <- !isFALSE(opts$cse)
  hoist <- !isFALSE(opts$hoist)
  dce <- !isFALSE(opts$dce)
  stmts <- model_statements(model)
  nodes <- unlist_nodes(lapply(stmts, function(st) parse_stmt(st$expr, st)))
  # linCmt() expands in place, into rxode2's compartments.
  lincmt <- NULL
  out <- list()
  for (n in nodes) {
    if (n$kind == "lincmt") {
      if (!is.null(lincmt)) abort_compile("linCmt()", n$st, "a second linCmt()")
      before <- assigned_names(out)
      lincmt <- lincmt_resolve(before, n$st)
      st <- n$st
      st$generated <- TRUE
      out <- c(out, lincmt_statements(lincmt, n$name, st))
    } else {
      out[[length(out) + 1L]] <- n
    }
  }
  nodes <- out
  lin_states <- if (is.null(lincmt)) character() else
    c(if (lincmt$depot) "depot", "central", if (lincmt$ncmt >= 2) "peripheral1",
      if (lincmt$ncmt >= 3) "peripheral2")
  states <- unique(c(lin_states, deriv_states(nodes)))
  if (!length(states)) {
    abort_compile("model", stmts[[1]], "a model with no compartment",
                  "d/dt(x) <- ... or linCmt()")
  }
  locals <- setdiff(assigned_names(nodes), states)
  clash <- intersect(assigned_names(nodes), states)
  if (length(clash)) {
    n <- Find(function(n) n$kind %in% c("assign") && n$name %in% clash, nodes)
    abort_compile(paste(clash[1], "<-"), n$st,
                  paste0("an assignment to the compartment `", clash[1], "`"),
                  paste0("d/dt(", clash[1], ") or ", clash[1], "(0)"))
  }
  params <- resolve_names(nodes, states, locals)
  class <- classify(nodes, states, !is.null(lincmt))
  nodes <- as_assignments(nodes)
  ny <- length(states)
  captures <- locals

  inputs <- c(
    stats::setNames(lapply(seq_along(states) - 1L, function(i) list(kind = "y", index = i)), states),
    stats::setNames(lapply(seq_along(params) - 1L, function(i) list(kind = "p", index = i)), params)
  )
  hoisted <- if (hoist) hoistable(nodes, states) else character()
  work <- stats::setNames(lapply(seq_along(hoisted) - 1L, function(i) list(kind = "w", index = i)), hoisted)
  # In rhs and post the hoisted locals are read from work, so their
  # statements are dropped there.
  drop_hoisted <- function(nodes) Filter(function(n) !(n$kind == "assign" && n$name %in% hoisted), nodes)

  pre_out <- unlist(lapply(pre_groups, function(g) paste0(".", g, ".", states)))
  defaults <- lapply(seq_along(pre_out), function(i) {
    g <- pre_groups[(i - 1L) %/% ny + 1L]
    list(kind = "assign", name = pre_out[i], expr = list(op = "num", value = pre_defaults[[g]]))
  })
  pre_nodes <- needed(c(defaults, nodes), pre_out, dce)
  # What pre's outputs depend on, whatever dce is: with dce off pre holds
  # statements that read states, whose values its outputs never use.
  pre_deps <- needed(c(defaults, nodes), pre_out, TRUE)
  pre_reads_state <- intersect(unique(unlist(lapply(pre_deps, node_reads))), states)
  if (length(pre_reads_state) || any(vapply(pre_deps, node_time, NA))) {
    st <- Find(function(n) any(startsWith(node_writes(n), ".")), Filter(function(n) n$kind != "assign" || !is.null(n$st), pre_deps))
    st <- if (is.null(st$st)) stmts[[1]] else st$st
    abort_compile("alag/f/rate/dur/x(0)", st,
                  "a lag, bioavailability, rate, duration or initial amount that depends on a compartment or on time",
                  "an expression of parameters and covariates only")
  }
  blocks <- list(
    pre = gen_block(pre_nodes, inputs, as.list(pre_out), fold = fold, cse = cse),
    ivl = if (length(hoisted)) gen_block(needed(nodes, hoisted, dce), inputs, as.list(hoisted),
                                         store_work = TRUE, fold = fold, cse = cse) else "",
    rhs = gen_block(needed(drop_hoisted(nodes), paste0(".d.", states), dce),
                    c(inputs, work), as.list(paste0(".d.", states)), fold = fold, cse = cse),
    post = gen_block(needed(drop_hoisted(nodes), captures, dce), c(inputs, work),
                     as.list(captures), fold = fold, cse = cse)
  )
  ir <- ir_asm(blocks, ny = ny, npar = length(params), nwork = length(hoisted),
               nout = c(pre = length(pre_out), rhs = ny, post = length(captures)))
  list(ir = ir, states = states, params = params, captures = captures,
       hoisted = hoisted, class = class, lincmt = lincmt,
       opts = list(fold = fold, cse = cse, hoist = hoist, dce = dce))
}
