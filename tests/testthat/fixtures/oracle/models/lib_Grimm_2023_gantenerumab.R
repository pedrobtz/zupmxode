cl <- exp(lcl + log(WT/5) * allo_cl)
vc <- exp(lvc + log(WT/5) * allo_v)
vp <- exp(lvp + log(WT/5) * allo_v)
q <- exp(lq + log(WT/5) * allo_cl)
fpla_cerebellum_i <- fpla_cerebellum * exp(bsv_fpla_cerebellum)
fpla_hippocampus_i <- fpla_hippocampus * exp(bsv_fpla_hippocampus)
fpla_striatum_i <- fpla_striatum * exp(bsv_fpla_striatum)
fpla_cortex_i <- fpla_cortex * exp(bsv_fpla_cortex)
fpla_choroid_plexus_i <- fpla_choroid_plexus * exp(bsv_fpla_choroid_plexus)
kel <- cl/vc
k12 <- q/vc
k21 <- q/vp
d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
d/dt(peripheral1) <- k12 * central - k21 * peripheral1
Cc <- central/vc * 1e+06
d/dt(Ccerebellum) <- kout_cerebellum * (kp_cerebellum * Cc - Ccerebellum)
Cbrain_cerebellum <- fpla_cerebellum_i * Cc + Ccerebellum
d/dt(Chippocampus) <- kout_hippocampus * (kp_hippocampus * Cc - Chippocampus)
Cbrain_hippocampus <- fpla_hippocampus_i * Cc + Chippocampus
d/dt(Cstriatum) <- kout_striatum * (kp_striatum * Cc - Cstriatum)
Cbrain_striatum <- fpla_striatum_i * Cc + Cstriatum
d/dt(Ccortex) <- kout_cortex * (kp_cortex * Cc - Ccortex)
Cbrain_cortex <- fpla_cortex_i * Cc + Ccortex
d/dt(Cchoroid_plexus) <- kout_choroid_plexus * (kp_choroid_plexus * Cc - Cchoroid_plexus)
Cbrain_choroid_plexus <- fpla_choroid_plexus_i * Cc + Cchoroid_plexus
d/dt(Ccsf) <- kout_csf * (kp_csf * Cc - Ccsf)
Cbrain_csf <- fpla_csf * Cc + Ccsf
