# Event tables (design §4.2, §6.2): NONMEM-style records, built by
# pmx_ode_events() or given as a data frame, checked, and scheduled by the
# C scheduler (src/core/sched.c) into one subject's timeline.

event_columns <- c("ID", "TIME", "EVID", "AMT", "CMT", "RATE", "DUR", "II",
                   "ADDL", "SS")

#' Build an event table
#'
#' Doses and observation times as NONMEM-style records: one dose record per
#' element of `amt` (recycled with `dose_time`, `cmt`, `rate`, `dur`, `ii`,
#' `addl` and `ss`), and one observation record per element of `time`.
#'
#' @param amt Dose amounts.
#' @param cmt The compartment each dose goes into, by name or number.
#' @param time Observation times.
#' @param dose_time When each dose is given.
#' @param ii,addl Interval and number of additional doses after each one.
#' @param rate,dur Infusion rate or duration; 0 for a bolus. `rate = -1`
#'   takes the rate from the model's `rate()`, `rate = -2` the duration
#'   from its `dur()`.
#' @param ss Steady state: 1 resets and doses to steady state, 2 adds the
#'   steady state to what is there.
#' @param id The subject.
#' @return A data frame of class `pmx_ode_events` with columns ID, TIME,
#'   EVID, AMT, CMT, RATE, DUR, II, ADDL and SS, sorted by time with doses
#'   before observations at the same time.
#' @examples
#' pmx_ode_events(amt = 100, ii = 24, addl = 6, cmt = "depot",
#'                time = seq(0, 168, by = 24))
#' @export
pmx_ode_events <- function(amt = double(), cmt = 1L, time = double(),
                           dose_time = 0, ii = 0, addl = 0L, rate = 0,
                           dur = 0, ss = 0L, id = 1L) {
  n <- if (length(amt)) max(length(amt), length(dose_time)) else 0L
  rec <- function(x) rep_len(x, n)
  doses <- data.frame(
    ID = rep_len(id, n), TIME = rec(as.double(dose_time)), EVID = rep_len(1L, n),
    AMT = rec(as.double(amt)), CMT = rec(cmt), RATE = rec(as.double(rate)),
    DUR = rec(as.double(dur)), II = rec(as.double(ii)),
    ADDL = rec(as.integer(addl)), SS = rec(as.integer(ss)),
    stringsAsFactors = FALSE
  )
  obs_cmt <- if (n) doses$CMT[1] else cmt[1]
  obs <- data.frame(
    ID = rep_len(id, length(time)), TIME = as.double(time),
    EVID = rep_len(0L, length(time)), AMT = 0, CMT = rep_len(obs_cmt, length(time)),
    RATE = 0, DUR = 0, II = 0, ADDL = 0L, SS = 0L, stringsAsFactors = FALSE
  )
  d <- rbind(doses, obs)
  d <- d[order(d$ID, d$TIME, d$EVID == 0L), , drop = FALSE]
  rownames(d) <- NULL
  structure(d, class = c("pmx_ode_events", "data.frame"))
}

abort_events <- function(row, rule, call = NULL) {
  pmx_abort("pmx_ode_error_events",
            paste0("Event table row ", row, ": ", rule, "."),
            row = row, rule = rule, call = call)
}

# The records of `events` as the scheduler reads them: NONMEM column names
# in any case, defaults for the ones left out, compartments by name mapped
# to the model's numbers. Covariates are the other numeric columns.
as_records <- function(events, model) {
  if (!is.data.frame(events)) {
    pmx_abort("pmx_ode_error_argument", "`events` must be a data frame.")
  }
  d <- as.data.frame(events, stringsAsFactors = FALSE)
  up <- toupper(names(d))
  hit <- up %in% event_columns
  names(d)[hit] <- up[hit]
  if (!"TIME" %in% names(d)) {
    pmx_abort("pmx_ode_error_events", "The event table has no TIME column.",
              row = NA_integer_, rule = "TIME is required")
  }
  n <- nrow(d)
  num <- function(col, default) {
    x <- if (col %in% names(d)) d[[col]] else rep(default, n)
    if (!is.numeric(x) && !is.logical(x)) {
      abort_events(1L, paste0(col, " must be numeric"))
    }
    x <- as.double(x)
    bad <- which(is.na(x))
    if (length(bad)) abort_events(bad[1], paste0(col, " is missing"))
    x
  }
  amt <- num("AMT", 0)
  evid <- if ("EVID" %in% names(d)) num("EVID", 0) else ifelse(amt > 0, 1, 0)
  cmt <- if ("CMT" %in% names(d)) d$CMT else rep(1L, n)
  if (is.factor(cmt)) cmt <- as.character(cmt)
  if (is.character(cmt)) {
    neg <- startsWith(cmt, "-")
    name <- sub("^-", "", cmt)
    idx <- match(name, model$states)
    bad <- which(is.na(idx))
    if (length(bad)) {
      abort_events(bad[1], paste0("CMT \"", cmt[bad[1]], "\" is not a compartment of the model (",
                                  paste(model$states, collapse = ", "), ")"))
    }
    cmt <- ifelse(neg, -idx, idx)
  } else {
    cmt <- num("CMT", 1)
  }
  covariates <- intersect(setdiff(names(d), c(event_columns, "ID")), model$params)
  list(
    time = num("TIME", 0), evid = as.integer(evid), amt = amt,
    cmt = as.integer(cmt), rate = num("RATE", 0), dur = num("DUR", 0),
    ii = num("II", 0), addl = as.integer(num("ADDL", 0)),
    ss = as.integer(num("SS", 0)), id = if ("ID" %in% names(d)) d$ID else rep(1L, n),
    covariates = d[covariates]
  )
}

# The parameter vector at each record: `params`, with each covariate column
# of the table taking its value at that record (last observation carried
# forward is the table's own: each record carries its value).
record_params <- function(model, rec, params) {
  # Missing means not given at all: a parameter given as NaN goes on, for
  # the scheduler or the solver to refuse where it matters.
  p <- params[model$params]
  missing <- model$params[!model$params %in% c(names(params), names(rec$covariates))]
  if (length(missing)) {
    pmx_abort("pmx_ode_error_argument",
              paste0("No value for the model's parameter", if (length(missing) > 1) "s",
                     " ", paste0("`", missing, "`", collapse = ", "), "."),
              call = NULL)
  }
  out <- matrix(rep(as.double(p), length(rec$time)), ncol = length(model$params),
                byrow = TRUE, dimnames = list(NULL, model$params))
  for (cv in names(rec$covariates)) out[, cv] <- as.double(rec$covariates[[cv]])
  out
}

# One subject's timeline, as a data frame: what the solvers will step
# through. For the tests (roadmap Stage 4).
sched_timeline <- function(model, events, params = double()) {
  rec <- as_records(events, model)
  if (length(unique(rec$id)) > 1L) {
    pmx_abort("pmx_ode_error_argument", "`events` must hold one subject.")
  }
  par <- record_params(model, rec, params)
  ny <- length(model$states)
  # The model's lag, bioavailability, rate and duration for each record's
  # compartment, from pre at that record's parameters.
  pre <- matrix(c(0, 1, 0, 0), nrow = length(rec$time), ncol = 4L, byrow = TRUE)
  for (i in which(rec$evid %in% c(1L, 4L))) {
    c <- rec$cmt[i]
    if (c < 1L || c > ny) next
    out <- model_eval(model, "pre", par = par[i, ], y = double(ny))$out
    g <- matrix(out, nrow = ny)
    pre[i, ] <- g[c, 1:4]
  }
  res <- .Call(zpx_schedule_r, rec$time, rec$evid, rec$amt, rec$cmt, rec$rate,
               rec$dur, rec$ii, rec$addl, rec$ss, pre[, 1], pre[, 2], pre[, 3],
               pre[, 4], ny)
  if (!is.null(res$error)) abort_events(res$row + 1L, res$error)
  kinds <- c("obs", "bolus", "inf_start", "inf_end", "reset", "ss", "off", "on", "mark")
  data.frame(time = res$time, kind = kinds[res$kind + 1L], cmt = res$cmt,
             amount = res$amount, rate = res$rate, dur = res$dur, ii = res$ii,
             ss = res$ss, row = res$rec + 1L, stringsAsFactors = FALSE)
}
