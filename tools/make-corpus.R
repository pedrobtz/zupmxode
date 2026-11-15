#!/usr/bin/env Rscript
# Writes the oracle fixtures every conformance test compares with (roadmap
# Stage 1, design §7): the paired corpus of tools/corpus/ solved by rxode2
# and mrgsolve, and nlmixr2lib's model library solved by rxode2.
#
#   Rscript tools/make-corpus.R            # from the package root
#
# Output, under tests/testthat/fixtures/oracle/ (regenerated, never edited):
#   models/<model>.R           the model text zupmxode will compile (§6.1)
#   models/<model>.params.csv  its parameters
#   events/<case>.csv          the event table, compartments by name
#   <engine>/<case>.csv        id, time and each output, at every
#                              observation record
#   manifest.csv               one row per case and engine: status, engine
#                              version, tolerances, date
#   disagreements.csv          cases where rxode2 and mrgsolve differ beyond
#                              the §6.3 tolerances; each needs a note in the
#                              hand-written disagreements.md
#
# A case is <model>__<scenario>, or lib_<name> for a library model, which
# has one scenario (paths stay under the tarball's 100-byte limit).
#
# Both engines solve at rtol 1e-10 and atol 1e-12, far inside the §6.3 tolerances zupmxode is compared at (rtol 1e-6,
# atol 1e-8), so a fixture's own solver error does not count against it.

suppressPackageStartupMessages({
  library(rxode2)
  library(mrgsolve)
  library(nlmixr2lib)
})

out <- file.path("tests", "testthat", "fixtures", "oracle")
stopifnot(file.exists("DESCRIPTION"))
source(file.path("tools", "corpus", "models.R"))
source(file.path("tools", "corpus", "scenarios.R"))

rtol <- 1e-10
atol <- 1e-12
cmp_rtol <- 1e-6
cmp_atol <- 1e-8
today <- format(Sys.Date())
versions <- c(rxode2 = as.character(packageVersion("rxode2")),
              mrgsolve = as.character(packageVersion("mrgsolve")),
              nlmixr2lib = as.character(packageVersion("nlmixr2lib")))

# Everything generated is replaced; disagreements.md, written by hand, stays.
unlink(file.path(out, c("models", "events", "rxode2", "mrgsolve")),
       recursive = TRUE)
unlink(file.path(out, c("manifest.csv", "disagreements.csv")))
for (d in c("models", "events", "rxode2", "mrgsolve")) {
  dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)
}

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, quote = TRUE)
}

# 17 significant digits, so a fixture read back is the double it was.
round_trip <- function(d) {
  for (n in names(d)) {
    if (is.double(d[[n]])) d[[n]] <- as.numeric(sprintf("%.17g", d[[n]]))
  }
  d
}

# ---------------------------------------------------------------------------
# Engines. Each takes a model, an event table with compartments by name and
# a parameter vector, and returns id, time and the outputs at every
# observation record, in record order.

# The rxode2 solve is the tests' own (tests/testthat/helper-oracle.R), so
# a fixture and the test that re-solves it cannot drift apart.
testthat_path <- file.path("tests", "testthat")
source(file.path(testthat_path, "helper-oracle.R"))
# Outside testthat, the fixtures are read from the source tree.
oracle_path <- function(...) file.path(testthat_path, "fixtures", "oracle", ...)
solve_rxode2 <- function(rx_text, events, params, outputs) {
  oracle_solve_rxode2(rx_text, events, params, outputs, rtol, atol)
}

mrg_cache <- new.env()
solve_mrgsolve <- function(name, mrg_text, events, params, outputs, states) {
  mod <- mrg_cache[[name]]
  if (is.null(mod)) {
    mod <- mrgsolve::mcode(paste0("corpus_", name), mrg_text, quiet = TRUE)
    assign(name, mod, envir = mrg_cache)
  }
  numbers <- match(events$CMT, names(mrgsolve::init(mod)))
  stopifnot(!anyNA(numbers))
  d <- data.frame(ID = events$ID, time = events$TIME, evid = events$EVID,
                  amt = events$AMT, cmt = numbers, rate = events$RATE,
                  ii = events$II, addl = events$ADDL, ss = events$SS)
  extra <- setdiff(names(events), c("ID", "TIME", "EVID", "AMT", "CMT",
                                    "RATE", "II", "ADDL", "SS"))
  for (e in extra) d[[e]] <- events[[e]]
  mod <- mrgsolve::update(mod, rtol = rtol, atol = atol, maxsteps = 1e6,
                          ss_rtol = rtol, ss_atol = atol)
  mod <- mrgsolve::param(mod, as.list(params))
  res <- as.data.frame(mrgsolve::mrgsim(mod, data = d, obsonly = TRUE,
                                        nocb = FALSE, digits = 17))
  data.frame(id = as.integer(res$ID), time = res$time,
             res[outputs], check.names = FALSE)
}

manifest <- list()
record <- function(case, engine, status, message = "") {
  manifest[[length(manifest) + 1L]] <<- data.frame(
    case = case, engine = engine, status = status,
    message = gsub("[\r\n]+", " ", message), version = versions[[engine]],
    rtol = rtol, atol = atol, date = today
  )
}

run <- function(case, engine, f) {
  res <- tryCatch(f(), error = function(e) e)
  if (inherits(res, "error")) {
    record(case, engine, "error", conditionMessage(res))
    return(NULL)
  }
  write_csv(round_trip(res), file.path(out, engine, paste0(case, ".csv")))
  record(case, engine, "ok")
  res
}

# The largest difference between two engines' tables beyond the §6.3
# tolerances, or NULL when they agree.
worst_gap <- function(a, b, outputs) {
  if (nrow(a) != nrow(b) || any(a$time != b$time)) {
    return(data.frame(output = "", time = NA, rxode2 = NA, mrgsolve = NA,
                      rel = NA, note = "different observation rows"))
  }
  worst <- NULL
  for (o in outputs) {
    gap <- abs(a[[o]] - b[[o]])
    scale <- pmax(abs(a[[o]]), abs(b[[o]]))
    bad <- which(gap > cmp_atol & gap > cmp_rtol * scale)
    if (length(bad)) {
      i <- bad[which.max(gap[bad] / pmax(scale[bad], cmp_atol))]
      worst <- rbind(worst, data.frame(
        output = o, time = a$time[i], rxode2 = a[[o]][i],
        mrgsolve = b[[o]][i], rel = gap[i] / max(scale[i], cmp_atol),
        note = sprintf("%d of %d rows", length(bad), nrow(a))
      ))
    }
  }
  worst
}

# ---------------------------------------------------------------------------
# The paired corpus.

disagreements <- list()
for (m in corpus_models) {
  writeLines(trimws(m$rx), file.path(out, "models", paste0(m$name, ".R")))
  write_csv(data.frame(name = names(m$params), value = unname(m$params)),
            file.path(out, "models", paste0(m$name, ".params.csv")))
  for (sname in names(corpus_scenarios)) {
    s <- corpus_scenarios[[sname]]
    if (!s$feature %in% m$features) next
    case <- paste0(m$name, "__", sname)
    events <- s$events(m)
    write_csv(events, file.path(out, "events", paste0(case, ".csv")))
    params <- m$params
    a <- run(case, "rxode2", function() {
      solve_rxode2(m$rx, events, params, m$outputs)
    })
    b <- run(case, "mrgsolve", function() {
      solve_mrgsolve(m$name, m$mrg, events, params, m$outputs, m$states)
    })
    if (!is.null(a) && !is.null(b)) {
      w <- worst_gap(a, b, m$outputs)
      if (!is.null(w)) disagreements[[case]] <- cbind(case = case, w)
    }
    cat(sprintf("%-40s %s\n", case,
                if (is.null(disagreements[[case]])) "" else "DISAGREE"))
  }
}

# ---------------------------------------------------------------------------
# nlmixr2lib: each library model without its residual-error lines (the
# §6.1 subset has none), at its initial thetas with every eta 0, dosed into
# its first compartment: 100 units every 24 time units, three doses.
# rxode2 is the only engine; covariates take the reference values below,
# and a model needing one not listed is recorded as skipped.

reference_covariates <- c(WT = 70, wt = 70, AGE = 40, age = 40, SEX = 0,
                          sex = 0, CRCL = 100, ALB = 4, BMI = 24)

for (name in nlmixr2lib::modeldb$name) {
  case <- paste0("lib_", name)
  ui <- tryCatch(rxode2::rxode2(nlmixr2lib::readModelDb(name)),
                 error = function(e) e)
  if (inherits(ui, "error")) {
    record(case, "rxode2", "skipped", conditionMessage(ui))
    next
  }
  lines <- ui$lstExpr
  keep <- !vapply(lines, function(e) is.call(e) &&
                    identical(e[[1]], as.name("~")), NA)
  rx_text <- paste(vapply(lines[keep], deparse1, ""), collapse = "\n")
  theta <- ui$theta
  etas <- ui$iniDf$name[!is.na(ui$iniDf$neta1) &
                          ui$iniDf$neta1 == ui$iniDf$neta2]
  params <- c(theta, stats::setNames(rep(0, length(etas)), etas))
  mod <- tryCatch(rxode2::rxode2(rx_text), error = function(e) e)
  if (inherits(mod, "error")) {
    record(case, "rxode2", "skipped", conditionMessage(mod))
    next
  }
  needed <- setdiff(rxode2::rxModelVars(mod)$params, names(params))
  if (length(setdiff(needed, names(reference_covariates)))) {
    record(case, "rxode2", "skipped", paste(
      "covariates without a reference value:",
      paste(setdiff(needed, names(reference_covariates)), collapse = ", ")
    ))
    next
  }
  params <- c(params, reference_covariates[needed])
  outputs <- intersect(c(ui$predDf$var, ui$predDf$cond), c(
    rxode2::rxModelVars(mod)$lhs, rxode2::rxModelVars(mod)$state
  ))
  if (!length(outputs)) {
    record(case, "rxode2", "skipped", "no prediction variable found")
    next
  }
  dose <- rxode2::rxModelVars(mod)$state[1]
  events <- ev_table(ev_row(0, amt = 100, cmt = dose, ii = 24, addl = 2),
                     ev_obs(c(0, 0.5, seq(1, 96, by = 1)), dose))
  writeLines(rx_text, file.path(out, "models", paste0("lib_", name, ".R")))
  write_csv(data.frame(name = names(params), value = unname(params)),
            file.path(out, "models", paste0("lib_", name, ".params.csv")))
  write_csv(events, file.path(out, "events", paste0(case, ".csv")))
  run(case, "rxode2", function() {
    solve_rxode2(rx_text, events, params, unique(outputs))
  })
  cat(sprintf("%-40s\n", case))
}

write_csv(do.call(rbind, manifest), file.path(out, "manifest.csv"))
d <- do.call(rbind, disagreements)
if (is.null(d)) {
  d <- data.frame(case = character(), output = character(), time = numeric(),
                  rxode2 = numeric(), mrgsolve = numeric(), rel = numeric(),
                  note = character())
}
write_csv(d, file.path(out, "disagreements.csv"))
m <- do.call(rbind, manifest)
cat(sprintf("\n%d cases: %d ok, %d error, %d skipped; %d disagreements\n",
            length(unique(m$case)), sum(m$status == "ok"),
            sum(m$status == "error"), sum(m$status == "skipped"), nrow(d)))
