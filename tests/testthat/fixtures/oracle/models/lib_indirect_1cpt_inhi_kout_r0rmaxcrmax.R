ka <- exp(lka)
vc <- exp(lvc)
cl <- exp(lcl)
r0 <- exp(lr0)
rmax <- exp(lrmax)
s1 <- exp(ls1)
imax <- exp(limax)
crmax <- exp(lcrmax)
fdepot <- exp(lfdepot)
kel <- cl/vc
imax <- (rmax - r0)/rmax
kin <- s1/imax
kout <- kin/r0
IC50 <- crmax * (r0 - (1 - imax) * rmax)/(rmax - r0)
Cc <- central/vc
d/dt(depot) <- -ka * depot
f(depot) <- fdepot
d/dt(central) <- ka * depot - (kel) * central
d/dt(effect) <- kin - kout * (1 - imax * Cc/(Cc + IC50)) * effect
