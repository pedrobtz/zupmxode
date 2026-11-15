vc <- exp(lvc + etalvc) * (WT/75)^e_wt_vc
ke <- exp(lke + etalke)
kcp <- exp(lkcp)
ka <- exp(lka + etalka)
MTT <- exp(lMTT + etamtt)
Vm <- exp(lVm + etalvm)
kpc <- kcp/Mpc
ktr <- (3 + 1)/MTT
d/dt(depot) <- -ktr * depot
d/dt(transit1) <- ktr * (depot - transit1)
d/dt(transit2) <- ktr * (transit1 - transit2)
d/dt(transit3) <- ktr * transit2 - ka * transit3
d/dt(central) <- ka * transit3 - ke * central - kcp * central + kpc * periph - central * (Vm/(Km + central/vc))
d/dt(periph) <- kcp * central - kpc * periph
f(depot) <- exp(lfdepot)
Cc <- central/vc
