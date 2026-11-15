ka <- exp(lka)
cl <- exp(lcl)
vc <- exp(lvc)
kel <- cl/vc
d/dt(depot) <- -ka * depot
d/dt(central) <- ka * depot - kel * central
Cc <- central/vc
