ka <- exp(lka)
vc <- exp(lvc)
cl <- exp(lcl)
IC50 <- exp(lIC50)
imax <- exp(limax)
kin <- exp(lkin)
kout <- exp(lkout)
fdepot <- exp(lfdepot)
kel <- cl/vc
d/dt(depot) <- -ka * depot
f(depot) <- fdepot
d/dt(central) <- ka * depot - (kel) * central
Cc <- central/vc
d/dt(effect) <- kin * (1 - imax * Cc/(Cc + IC50)) - kout * effect
