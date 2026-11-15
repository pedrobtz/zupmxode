d/dt(depot) <- -ka * depot
d/dt(central) <- ka * depot - (cl + q) / v * central + q / vp * periph
d/dt(periph) <- q / v * central - q / vp * periph
cp <- central / v
