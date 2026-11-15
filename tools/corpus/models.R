# The paired conformance corpus: each model written twice, in the rxode2
# subset of design §6.1 and in mrgsolve's language, with one parameter set.
# tools/make-corpus.R solves every model with both engines on the
# scenarios of tools/corpus/scenarios.R. Models are ports of mrgsolve's
# modlib() (pk1, pk2, pk3, irm1, emax) and of the event features §6.2
# names; nlmixr2lib's library is added by make-corpus.R itself.
#
# Fields: rx and mrg are model text; params is shared by both; states are
# the compartments in the order both engines number them; dose is the
# compartment doses go to by default; outputs are the captured values
# compared; class is the solver class design §4.1 expects; features lists
# the scenario features (tools/corpus/scenarios.R) the model supports.

corpus_model <- function(name, rx, mrg, params, states, dose, outputs = "cp",
                         class, features = character(), covariates = NULL) {
  list(name = name, rx = rx, mrg = mrg, params = params, states = states,
       dose = dose, outputs = outputs, class = class, features = features,
       covariates = covariates)
}

mrg_param <- function(params) {
  paste0("$PARAM ", paste(names(params), "=", params, collapse = ", "))
}

pk1 <- c(ka = 0.8, cl = 2, v = 20)
pk2 <- c(ka = 0.8, cl = 2, v = 20, q = 4, vp = 60)
pk3 <- c(cl = 2, v = 20, q = 4, vp = 60, q2 = 1, vp2 = 100)

# linCmt() finds its parameters among the variables a model assigns, so the
# closed-form models take t-prefixed parameters and assign them.
tpar <- function(p) stats::setNames(p, paste0("t", names(p)))
assign_lines <- function(p) paste0(names(p), " <- t", names(p), collapse = "\n")

corpus_models <- list(
  corpus_model(
    "pk1_oral", class = "linear", params = pk1,
    states = c("depot", "central"), dose = "depot",
    features = c("bolus", "addl", "ss", "reset"),
    rx = "
d/dt(depot) <- -ka * depot
d/dt(central) <- ka * depot - cl / v * central
cp <- central / v
",
    mrg = paste(mrg_param(pk1), "
$CMT depot central
$ODE
dxdt_depot = -ka * depot;
dxdt_central = ka * depot - cl / v * central;
$TABLE capture cp = central / v;
")
  ),
  corpus_model(
    "pk1_oral_lin", class = "closed-form", params = tpar(pk1),
    states = c("depot", "central"), dose = "depot",
    features = c("bolus", "addl", "ss", "reset"),
    rx = paste0(assign_lines(pk1), "\ncp <- linCmt()\n"),
    mrg = paste(mrg_param(tpar(pk1)), "
$CMT depot central
$PKMODEL ncmt = 1, depot = TRUE
$MAIN
double CL = tcl;
double V = tv;
double KA = tka;
$TABLE capture cp = central / V;
")
  ),
  corpus_model(
    "pk1_iv", class = "linear", params = pk1[c("cl", "v")],
    states = "central", dose = "central",
    features = c("bolus", "addl", "ss", "reset", "infusion"),
    rx = "
d/dt(central) <- -cl / v * central
cp <- central / v
",
    mrg = paste(mrg_param(pk1[c("cl", "v")]), "
$CMT central
$ODE dxdt_central = -cl / v * central;
$TABLE capture cp = central / v;
")
  ),
  corpus_model(
    "pk1_iv_ratedur", class = "linear",
    params = c(pk1[c("cl", "v")], r0 = 50, d0 = 3),
    states = "central", dose = "central",
    features = c("model_rate", "model_dur"),
    rx = "
d/dt(central) <- -cl / v * central
rate(central) <- r0
dur(central) <- d0
cp <- central / v
",
    mrg = paste(mrg_param(c(pk1[c("cl", "v")], r0 = 50, d0 = 3)), "
$CMT central
$MAIN
R_central = r0;
D_central = d0;
$ODE dxdt_central = -cl / v * central;
$TABLE capture cp = central / v;
")
  ),
  corpus_model(
    "pk1_oral_lagf", class = "linear",
    params = c(pk1, lag = 1.5, fdepot = 0.8),
    states = c("depot", "central"), dose = "depot",
    features = c("bolus", "addl", "lag"),
    rx = "
alag(depot) <- lag
f(depot) <- fdepot
d/dt(depot) <- -ka * depot
d/dt(central) <- ka * depot - cl / v * central
cp <- central / v
",
    mrg = paste(mrg_param(c(pk1, lag = 1.5, fdepot = 0.8)), "
$CMT depot central
$MAIN
ALAG_depot = lag;
F_depot = fdepot;
$ODE
dxdt_depot = -ka * depot;
dxdt_central = ka * depot - cl / v * central;
$TABLE capture cp = central / v;
")
  ),
  corpus_model(
    "pk1_iv_lag", class = "linear", params = c(pk1[c("cl", "v")], lag = 1),
    states = "central", dose = "central",
    features = c("lag_infusion"),
    rx = "
alag(central) <- lag
d/dt(central) <- -cl / v * central
cp <- central / v
",
    mrg = paste(mrg_param(c(pk1[c("cl", "v")], lag = 1)), "
$CMT central
$MAIN ALAG_central = lag;
$ODE dxdt_central = -cl / v * central;
$TABLE capture cp = central / v;
")
  ),
  corpus_model(
    "pk1_iv_wt", class = "linear", params = c(pk1[c("cl", "v")], wt = 70),
    states = "central", dose = "central", covariates = "wt",
    features = c("covariate"),
    rx = "
cli <- cl * (wt / 70)^0.75
d/dt(central) <- -cli / v * central
cp <- central / v
",
    mrg = paste(mrg_param(c(pk1[c("cl", "v")], wt = 70)), "
$CMT central
$ODE
double cli = cl * pow(wt / 70.0, 0.75);
dxdt_central = -cli / v * central;
$TABLE capture cp = central / v;
")
  ),
  corpus_model(
    "pk2_oral", class = "linear", params = pk2,
    states = c("depot", "central", "periph"), dose = "depot",
    features = c("bolus", "addl", "ss", "reset"),
    rx = "
d/dt(depot) <- -ka * depot
d/dt(central) <- ka * depot - (cl + q) / v * central + q / vp * periph
d/dt(periph) <- q / v * central - q / vp * periph
cp <- central / v
",
    mrg = paste(mrg_param(pk2), "
$CMT depot central periph
$ODE
dxdt_depot = -ka * depot;
dxdt_central = ka * depot - (cl + q) / v * central + q / vp * periph;
dxdt_periph = q / v * central - q / vp * periph;
$TABLE capture cp = central / v;
")
  ),
  corpus_model(
    "pk2_iv_lin", class = "closed-form", params = tpar(pk2[-1]),
    states = c("central", "periph"), dose = "central",
    features = c("bolus", "addl", "ss", "reset", "infusion"),
    rx = paste0(assign_lines(pk2[-1]), "\ncp <- linCmt()\n"),
    mrg = paste(mrg_param(tpar(pk2[-1])), "
$CMT central periph
$PKMODEL ncmt = 2, depot = FALSE
$MAIN
double CL = tcl;
double V1 = tv;
double Q = tq;
double V2 = tvp;
$TABLE capture cp = central / V1;
")
  ),
  corpus_model(
    "pk3_iv_lin", class = "closed-form", params = tpar(pk3),
    states = c("central", "periph1", "periph2"), dose = "central",
    features = c("bolus", "addl", "ss", "infusion"),
    # mrgsolve has no three-compartment closed form; its side is the ODE.
    rx = paste0(assign_lines(pk3), "\ncp <- linCmt()\n"),
    mrg = paste(mrg_param(tpar(pk3)), "
$MAIN
double cl = tcl;
double v = tv;
double q = tq;
double vp = tvp;
double q2 = tq2;
double vp2 = tvp2;
$CMT central periph1 periph2
$ODE
dxdt_central = -(cl + q + q2) / v * central + q / vp * periph1 +
  q2 / vp2 * periph2;
dxdt_periph1 = q / v * central - q / vp * periph1;
dxdt_periph2 = q2 / v * central - q2 / vp2 * periph2;
$TABLE capture cp = central / v;
")
  ),
  corpus_model(
    "pk1_mm", class = "general", params = c(vmax = 10, km = 2, v = 20),
    states = "central", dose = "central",
    features = c("bolus", "addl", "ss", "infusion"),
    rx = "
cp <- central / v
d/dt(central) <- -vmax * cp / (km + cp)
",
    mrg = paste(mrg_param(c(vmax = 10, km = 2, v = 20)), "
$CMT central
$ODE
double c = central / v;
dxdt_central = -vmax * c / (km + c);
$TABLE capture cp = central / v;
")
  ),
  corpus_model(
    "irm1", class = "general",
    params = c(pk1, kin = 10, kout = 0.5, imax = 0.9, ic50 = 2),
    states = c("depot", "central", "resp"), dose = "depot",
    outputs = c("cp", "resp"), features = c("bolus", "addl"),
    rx = "
resp(0) <- kin / kout
cp <- central / v
d/dt(depot) <- -ka * depot
d/dt(central) <- ka * depot - cl / v * central
d/dt(resp) <- kin * (1 - imax * cp / (ic50 + cp)) - kout * resp
",
    mrg = paste(mrg_param(c(pk1, kin = 10, kout = 0.5, imax = 0.9,
                            ic50 = 2)), "
$CMT depot central resp
$MAIN resp_0 = kin / kout;
$ODE
double c = central / v;
dxdt_depot = -ka * depot;
dxdt_central = ka * depot - cl / v * central;
dxdt_resp = kin * (1 - imax * c / (ic50 + c)) - kout * resp;
$TABLE capture cp = central / v;
")
  ),
  corpus_model(
    "emax_effect", class = "linear",
    params = c(pk1, ke0 = 0.3, e0 = 5, emax = 40, ec50 = 1.5),
    states = c("depot", "central", "ce"), dose = "depot",
    outputs = c("cp", "eff"), features = c("bolus", "addl"),
    rx = "
d/dt(depot) <- -ka * depot
d/dt(central) <- ka * depot - cl / v * central
cp <- central / v
d/dt(ce) <- ke0 * (cp - ce)
eff <- e0 + emax * ce / (ec50 + ce)
",
    mrg = paste(mrg_param(c(pk1, ke0 = 0.3, e0 = 5, emax = 40,
                            ec50 = 1.5)), "
$CMT depot central ce
$ODE
double c = central / v;
dxdt_depot = -ka * depot;
dxdt_central = ka * depot - cl / v * central;
dxdt_ce = ke0 * (c - ce);
$TABLE
capture cp = central / v;
capture eff = e0 + emax * ce / (ec50 + ce);
")
  )
)
names(corpus_models) <- vapply(corpus_models, `[[`, "", "name")
