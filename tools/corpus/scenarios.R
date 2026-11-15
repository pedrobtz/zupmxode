# Event scenarios for the paired corpus: one per event rule of design §6.2
# that the corpus models can show. Each builds a NONMEM-style table, with
# compartments by name, for one model; tools/make-corpus.R maps the names
# to each engine's numbers. A scenario applies to every model whose
# features include its `feature`.
#
# Rows at the same time are ordered doses first, then observations, so an
# observation at a dose time is read after the dose, as NONMEM reads a data
# set sorted that way.

ev_row <- function(time, evid = 1, amt = 0, cmt, rate = 0, ii = 0, addl = 0,
                   ss = 0) {
  data.frame(ID = 1L, TIME = time, EVID = evid, AMT = amt, CMT = cmt,
             RATE = rate, II = ii, ADDL = addl, SS = ss)
}

ev_obs <- function(times, cmt) {
  ev_row(times, evid = 0, cmt = cmt)
}

ev_table <- function(...) {
  d <- do.call(rbind, list(...))
  d <- d[order(d$TIME, d$EVID == 0), , drop = FALSE]
  rownames(d) <- NULL
  d
}

obs_grid <- c(0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 36, 48)

corpus_scenarios <- list(
  bolus_single = list(feature = "bolus", events = function(m) {
    ev_table(ev_row(0, amt = 100, cmt = m$dose), ev_obs(obs_grid, m$dose))
  }),
  bolus_addl = list(feature = "addl", events = function(m) {
    ev_table(ev_row(0, amt = 100, cmt = m$dose, ii = 12, addl = 5),
             ev_obs(c(seq(0, 96, by = 2), 60.5, 61), m$dose))
  }),
  ss1_bolus = list(feature = "ss", events = function(m) {
    ev_table(ev_row(0, amt = 100, cmt = m$dose, ii = 12, ss = 1),
             ev_obs(c(0, 0.5, 1, 2, 4, 6, 8, 11.99, 12, 18, 24), m$dose))
  }),
  ss2_add = list(feature = "ss", events = function(m) {
    ev_table(ev_row(0, amt = 100, cmt = m$dose, ii = 24, ss = 1),
             ev_row(12, amt = 50, cmt = m$dose, ii = 24, ss = 2),
             ev_obs(c(0, 2, 6, 11.99, 12, 14, 18, 23.99, 24, 36), m$dose))
  }),
  reset_evid3 = list(feature = "reset", events = function(m) {
    ev_table(ev_row(0, amt = 100, cmt = m$dose),
             ev_row(24, evid = 3, cmt = m$dose),
             ev_row(30, amt = 100, cmt = m$dose),
             ev_obs(c(0, 2, 12, 23.99, 24, 26, 29.99, 30, 32, 48), m$dose))
  }),
  reset_dose_evid4 = list(feature = "reset", events = function(m) {
    ev_table(ev_row(0, amt = 100, cmt = m$dose),
             ev_row(24, evid = 4, amt = 50, cmt = m$dose),
             ev_obs(c(0, 2, 12, 23.99, 24, 26, 36, 48), m$dose))
  }),
  infusion_rate = list(feature = "infusion", events = function(m) {
    ev_table(ev_row(0, amt = 100, cmt = m$dose, rate = 25),
             ev_obs(c(obs_grid, 3.99, 4.01), m$dose))
  }),
  ss1_infusion = list(feature = "infusion", events = function(m) {
    ev_table(ev_row(0, amt = 100, cmt = m$dose, rate = 20, ii = 24, ss = 1),
             ev_obs(c(0, 1, 4.99, 5, 8, 12, 23.99, 24), m$dose))
  }),
  rate_modelled = list(feature = "model_rate", events = function(m) {
    ev_table(ev_row(0, amt = 100, cmt = m$dose, rate = -1),
             ev_obs(c(obs_grid, 1.99, 2.01), m$dose))
  }),
  dur_modelled = list(feature = "model_dur", events = function(m) {
    ev_table(ev_row(0, amt = 100, cmt = m$dose, rate = -2),
             ev_obs(c(obs_grid, 2.99, 3.01), m$dose))
  }),
  lag_bolus = list(feature = "lag", events = function(m) {
    ev_table(ev_row(0, amt = 100, cmt = m$dose, ii = 12, addl = 1),
             ev_obs(c(obs_grid, 1.49, 1.51, 13.49, 13.51), m$dose))
  }),
  # A lagged infusion whose end falls on the next infusion's start: the
  # first runs 1 to 3 after its lag, the second starts at 3 after its own.
  lag_infusion_straddle = list(feature = "lag_infusion", events = function(m) {
    ev_table(ev_row(0, amt = 100, cmt = m$dose, rate = 50),
             ev_row(2, amt = 100, cmt = m$dose, rate = 50),
             ev_obs(c(0, 0.5, 0.99, 1.01, 2, 2.99, 3, 3.01, 4, 4.99, 5,
                      5.01, 8, 12), m$dose))
  }),
  # A covariate that changes on an EVID 2 record: carried forward from
  # that record, not back from it (last observation carried forward).
  covariate_locf = list(feature = "covariate", events = function(m) {
    d <- ev_table(ev_row(0, amt = 100, cmt = m$dose),
                  ev_row(12, evid = 2, cmt = m$dose),
                  ev_obs(c(obs_grid, 11.99, 12.01), m$dose))
    d$wt <- ifelse(d$TIME < 12, 70, 50)
    d
  })
)
