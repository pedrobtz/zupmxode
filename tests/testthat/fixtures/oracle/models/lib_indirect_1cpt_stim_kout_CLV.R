ka <- exp(lka)
vc <- exp(lvc)
cl <- exp(lcl)
EC50 <- exp(lEC50)
Emax <- exp(lEmax)
kin <- exp(lkin)
kout <- exp(lkout)
fdepot <- exp(lfdepot)
kel <- cl/vc
d/dt(depot) <- -ka * depot
f(depot) <- fdepot
d/dt(central) <- ka * depot - (kel) * central
Cc <- central/vc
d/dt(effect) <- kin - kout * (1 + Emax * Cc/(Cc + EC50)) * effect
