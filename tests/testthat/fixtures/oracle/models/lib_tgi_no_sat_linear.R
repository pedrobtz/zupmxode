ts0 <- exp(lts0)
ka <- exp(lka)
cl <- exp(lcl)
vc <- exp(lvc)
kgl <- exp(lkgl)
kel <- cl/vc
tumorSize(0) <- ts0
d/dt(depot) <- -ka * depot
d/dt(central) <- ka * depot - kel * central
d/dt(tumorSize) <- kgl
Cc <- central/vc
