ka <- exp(lka)
vc <- exp(lvc)
cl <- exp(lcl)
IC50 <- exp(lIC50)
imax <- exp(limax)
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
d/dt(precursor) <- kpin - (kin + kpout) * (1 - imax * Cc/(Cc + IC50)) * precursor
d/dt(effect) <- kin * (1 - imax * Cc/(Cc + IC50)) * precursor - kout * effect
