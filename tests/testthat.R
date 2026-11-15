# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional testing?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(zupmxode)

# Shuffled, so a test that leans on another's leftovers fails: the tests in
# each file run in a random order. The seed is printed; ZUPMXODE_TEST_SEED
# replays an order.
seed <- suppressWarnings(as.integer(Sys.getenv("ZUPMXODE_TEST_SEED")))
if (is.na(seed)) seed <- as.integer(Sys.time()) %% 100000L
cat("test order seed:", seed, "(ZUPMXODE_TEST_SEED to replay)\n")
set.seed(seed)
test_check("zupmxode", shuffle = TRUE)
