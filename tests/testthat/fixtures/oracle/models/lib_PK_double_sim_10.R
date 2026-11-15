ka1 <- exp(lka1)
k02 <- exp(tk02)
cl <- exp(lcl)
vc <- exp(lvc)
fdepot1 <- expit(lgfdepot1, 0, 1)
alag <- exp(lalag)
kel <- cl/vc
d/dt(depot1) <- -ka1 * depot1
f(depot1) <- fdepot1
d/dt(depot2) <- -k02
lag(depot2) <- alag
f(depot2) <- 1 - fdepot1
d/dt(central) <- ka1 * depot1 + k02 - kel * central
Cc <- central/vc
