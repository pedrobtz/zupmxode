clts <- exp(lcl)
cltmax <- exp(lcltmax)
clgamma <- exp(lclgamma)
clt50 <- exp(lclt50)
vc <- exp(lvc)
vp <- exp(lvp)
q <- exp(lq)
cl <- clts * exp(cltmax * time^clgamma/(clt50^clgamma + time^clgamma))
kel <- cl/vc
k12 <- q/vc
k21 <- q/vp
d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
d/dt(peripheral1) <- k12 * central - k21 * peripheral1
Cc <- central/vc
