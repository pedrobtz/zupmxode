ts0 <- exp(lts0)
tsmax <- exp(ltsmax)
ka <- exp(lka)
cl <- exp(lcl)
vc <- exp(lvc)
kg <- exp(lkg)
kd <- exp(lkd)
gamma <- exp(lgamma)
kel <- cl/vc
tumorSize(0) <- ts0
d/dt(depot) <- -ka * depot
d/dt(central) <- ka * depot - kel * central
d/dt(tumorSize) <- kg * tumorSize^(gamma) - kd * tumorSize
Cc <- central/vc
