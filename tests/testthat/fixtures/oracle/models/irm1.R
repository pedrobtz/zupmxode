resp(0) <- kin / kout
cp <- central / v
d/dt(depot) <- -ka * depot
d/dt(central) <- ka * depot - cl / v * central
d/dt(resp) <- kin * (1 - imax * cp / (ic50 + cp)) - kout * resp
