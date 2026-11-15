ts0 <- exp(lts0)
tsmax <- exp(ltsmax)
ka <- exp(lka)
cl <- exp(lcl)
vc <- exp(lvc)
kgl <- exp(lkgl)
alpha <- exp(lalpha)
beta <- exp(lbeta)
kel <- cl/vc
tumorSize(0) <- ts0
d/dt(depot) <- -ka * depot
d/dt(central) <- ka * depot - kel * central
d/dt(tumorSize) <- (alpha - beta * log(tumorSize)) * tumorSize
Cc <- central/vc
