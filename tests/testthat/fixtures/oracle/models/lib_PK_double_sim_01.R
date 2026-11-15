k01 <- exp(tk01)
ka2 <- exp(lka2)
cl <- exp(lcl)
vc <- exp(lvc)
fdepot1 <- expit(lgfdepot1, 0, 1)
alag <- exp(lalag)
kel <- cl/vc
d/dt(depot1) <- -k01
f(depot1) <- fdepot1
d/dt(depot2) <- -ka2 * depot2
lag(depot2) <- alag
f(depot2) <- 1 - fdepot1
d/dt(central) <- k01 + ka2 * depot2 - kel * central
Cc <- central/vc
