cl <- exp(lcl)
vc <- exp(lvc)
vp <- exp(lvp)
q <- exp(lq)
kel <- cl/vc
k12 <- q/vc
k21 <- q/vp
d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
d/dt(peripheral1) <- k12 * central - k21 * peripheral1
Cc <- central/vc
