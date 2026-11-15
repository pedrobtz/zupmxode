d/dt(depot) <- -ka * depot
d/dt(central) <- ka * depot - cl / v * central
cp <- central / v
d/dt(ce) <- ke0 * (cp - ce)
eff <- e0 + emax * ce / (ec50 + ce)
