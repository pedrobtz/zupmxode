fdepot <- exp(lfdepot + etafdepot)
ka <- exp(lka + etaka)
wtnorm <- log(WT/70)
cl <- exp(lcl + allocl * wtnorm + etacl)
q <- exp(lq + allocl * wtnorm + etaq)
v <- exp(lv + allov * wtnorm + etav)
vp <- exp(lvp + allov * wtnorm + etavp)
Cc <- linCmt()
f(depot) <- fdepot
