# The classifier (design §4.1) against labels written by reading each model;
# the labels and their reasons are in helper-compile.R.


test_that("each corpus model's class is its hand label", {
  texts <- corpus_model_texts()
  for (name in names(hand_labels)) {
    expect_identical(pmx_ode_model(texts[[name]])$class, hand_labels[[name]],
                     label = name)
  }
  # Every corpus model the compiler accepts has a label.
  accepted <- setdiff(names(texts), "lib_tgi_no_sat_expo")
  expect_setequal(names(hand_labels), accepted)
})

test_that("each rule of the classifier", {
  cls <- function(text) pmx_ode_model(text)$class
  expect_identical(cls("d/dt(x) <- -k * x"), "linear")
  expect_identical(cls("d/dt(x) <- -k * x + r"), "linear")
  expect_identical(cls("d/dt(x) <- -x / v"), "linear")
  expect_identical(cls("d/dt(x) <- -x * x"), "general")
  expect_identical(cls("d/dt(x) <- -v / x"), "general")
  expect_identical(cls("d/dt(x) <- -exp(x)"), "general")
  expect_identical(cls("d/dt(x) <- -k * t * x"), "general")
  expect_identical(cls("c <- x / v\nd/dt(x) <- -k * c"), "linear")
  expect_identical(cls("c <- x / v\nd/dt(x) <- -c * c"), "general")
  expect_identical(cls("if (k > 1) a <- 1 else a <- 2\nd/dt(x) <- -a * x"), "linear")
  expect_identical(cls("if (x > 1) a <- 1 else a <- 2\nd/dt(x) <- -a * x"), "general")
  expect_identical(cls("cl <- 1\nv <- 2\ncp <- linCmt()"), "closed-form")
  expect_identical(cls("cl <- 1\nv <- 2\ncp <- linCmt()\nd/dt(e) <- k * (cp - e)"), "linear")
  expect_identical(cls("cl <- 1\nv <- 2\ncp <- linCmt()\nd/dt(e) <- k * cp * e"), "general")
})
