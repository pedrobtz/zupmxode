ka <- exp(lka)
vc <- exp(lvc)
cl <- exp(lcl)
EC50 <- exp(lEC50)
emax <- exp(lemax)
kin <- exp(lkin)
kout <- exp(lkout)
kpin <- exp(lkpin)
kpout <- exp(lkpout)
fdepot <- exp(lfdepot)
kel <- cl/vc
Cc <- central/vc
d/dt(depot) <- -ka * depot
f(depot) <- fdepot
d/dt(central) <- ka * depot - kel * central
d/dt(precursor) <- kpin - (kin + kpout) * (1 + emax * Cc/(Cc + EC50)) * precursor
d/dt(effect) <- kin * (1 + emax * Cc/(Cc + EC50)) * precursor - kout * effect
