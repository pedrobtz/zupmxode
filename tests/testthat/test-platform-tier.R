# The README says which tier each r-actions leg runs (roadmap Stage 11):
# native code and MIR's interpreter, both, on every leg of R-CMD-check.
# This test holds the README to it; a leg where either is missing fails
# here, and the README is corrected with it.

test_that("both tiers run on every R-CMD-check leg", {
  skip_if_not(identical(Sys.getenv("R_ACTIONS_WORKFLOW"), "r-cmd-check"),
              "the claim is about the CI legs")
  expect_identical(available_tiers(), c("native", "interp"))
  expect_identical(zumir::zumir_tier(), "native")
})
