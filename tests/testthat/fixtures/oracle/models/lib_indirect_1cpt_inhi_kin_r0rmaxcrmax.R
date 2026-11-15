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
imax <- (r0 - rmax)/r0
kin <- -s1/imax
kout <- kin/r0
IC50 <- crmax * (rmax - (1 - imax) * r0)/(r0 - rmax)
d/dt(depot) <- -ka * depot
f(depot) <- fdepot
d/dt(central) <- ka * depot - (kel) * central
Cc <- central/vc
d/dt(effect) <- kin * (1 - imax * Cc/(Cc + IC50)) - kout * effect
