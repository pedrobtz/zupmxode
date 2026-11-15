test_that("the R glue reaches the C core and reports its ABI", {
  expect_identical(.Call(zupmxode_core_abi_r), 1L)
})
