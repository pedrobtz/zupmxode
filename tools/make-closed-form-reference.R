#!/usr/bin/env Rscript
# Writes tests/testthat/fixtures/closed-form-reference.csv: amounts after
# 1,000 doses, by superposition at 200 bits (Rmpfr), the oracle for the
# closed-form solver's state propagation over a long dosing history
# (design D21). Superposition is exact here: the parameters never change.
#
#   Rscript tools/make-closed-form-reference.R      # needs Rmpfr
suppressPackageStartupMessages(library(Rmpfr))
stopifnot(file.exists("DESCRIPTION"))
M <- function(x) mpfr(x, 200)
bits <- function(x) {
  b <- matrix(writeBin(as.double(x), raw(), endian = "big"), nrow = 8)
  apply(b, 2, function(col) paste(format(col), collapse = ""))
}
doses <- (0:999) * 1           # 100 every hour, 1,000 doses
times <- c(0.5, 499.5, 999.5, 1000.5, 1100)
rows <- list()

# One compartment with a depot: central(t) = D ka / (ka - k) (e^-kt - e^-ka t).
ka <- M(0.8); k <- M(2) / M(20)
for (t in times) {
  dt <- M(t) - M(doses[doses <= t])
  c1 <- sum(100 * ka / (ka - k) * (exp(-k * dt) - exp(-ka * dt)))
  rows[[length(rows) + 1L]] <- data.frame(model = "pk1_oral", time = t, state = "central",
                                          ref = bits(asNumeric(c1)))
}

# Two compartments, IV bolus: central(t) = D (a e^-alpha t + b e^-beta t).
k10 <- M(2) / M(20); k12 <- M(4) / M(20); k21 <- M(4) / M(60)
s <- k10 + k12 + k21
alpha <- (s + sqrt(s^2 - 4 * k10 * k21)) / 2
beta <- (s - sqrt(s^2 - 4 * k10 * k21)) / 2
a <- (alpha - k21) / (alpha - beta)
b <- (k21 - beta) / (alpha - beta)
for (t in times) {
  dt <- M(t) - M(doses[doses <= t])
  c2 <- sum(100 * (a * exp(-alpha * dt) + b * exp(-beta * dt)))
  rows[[length(rows) + 1L]] <- data.frame(model = "pk2_iv", time = t, state = "central",
                                          ref = bits(asNumeric(c2)))
}
out <- file.path("tests", "testthat", "fixtures", "closed-form-reference.csv")
utils::write.csv(do.call(rbind, rows), out, row.names = FALSE)
cat("wrote", out, "\n")
