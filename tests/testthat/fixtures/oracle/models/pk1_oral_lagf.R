alag(depot) <- lag
f(depot) <- fdepot
d/dt(depot) <- -ka * depot
d/dt(central) <- ka * depot - cl / v * central
cp <- central / v
