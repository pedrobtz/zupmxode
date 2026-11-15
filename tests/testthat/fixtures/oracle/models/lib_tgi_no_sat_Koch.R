ts0 <- exp(lts0)
ka <- exp(lka)
cl <- exp(lcl)
vc <- exp(lvc)
kge <- exp(lkge)
kgl <- exp(lkgl)
kel <- cl/vc
tumorSize(0) <- ts0
d/dt(depot) <- -ka * depot
d/dt(central) <- ka * depot - kel * central
d/dt(tumorSize) <- (2 * kge * kgl * tumorSize)/(kgl + 2 * kge * tumorSize)
Cc <- central/vc
