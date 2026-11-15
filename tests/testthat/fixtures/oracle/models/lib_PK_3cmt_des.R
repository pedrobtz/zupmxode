ka <- exp(lka)
cl <- exp(lcl)
vc <- exp(lvc)
vp <- exp(lvp)
vp2 <- exp(lvp2)
q <- exp(lq)
q2 <- exp(lq2)
kel <- cl/vc
k12 <- q/vc
k21 <- q/vp
k13 <- q2/vc
k31 <- q2/vp2
d/dt(depot) <- -ka * depot
d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
d/dt(peripheral1) <- k12 * central - k21 * peripheral1
d/dt(peripheral2) <- k13 * central - k31 * peripheral2
Cc <- central/vc
