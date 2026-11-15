# Compile errors (design §8): a construct outside the subset is an error of
# class pmx_ode_error_compile naming the construct, its position and the
# nearest supported alternative.


test_that("rxode2's features outside the subset name themselves", {
  base <- "d/dt(x) <- -k * x\n"
  expect_compile_error(paste0(base, "y <- ifelse(x > 1, 1, 0)"), "ifelse", "if (...)")
  expect_compile_error(paste0(base, "y <- x\ny ~ prop(0.1)"), "~", "pmx_ode_sim(sigma = )")
  expect_compile_error(paste0(base, "for (i in 1:3) y <- i"), "for")
  expect_compile_error(paste0(base, "while (x > 1) y <- 1"), "while")
  expect_compile_error(paste0(base, "y <- transit(3, 1)"), "transit", "d/dt()")
  expect_compile_error(paste0(base, "mtime(t1) <- 2"), "mtime(t1) <-")
  expect_compile_error(paste0(base, "y <- sum(x, 1)"), "sum", "+")
  expect_compile_error(paste0(base, "y <- log10(x)"), "log10", "log(x) / log(10)")
  expect_compile_error(paste0(base, "y <- rnorm(1)"), "rnorm", "pmx_ode_sim()")
  expect_compile_error(paste0(base, "y <<- 1"), "<<-", "<-")
  expect_compile_error(paste0(base, "y <- x[1]"), "[")
})

test_that("a misspelt function suggests the table's name", {
  expect_compile_error("d/dt(x) <- -expp(k) * x", "expp", "exp()")
  expect_compile_error("d/dt(x) <- -exp(k, 2) * x", "exp()", "exp(x)")
})

test_that("misuse of the language is refused where it is written", {
  err <- expect_compile_error("d/dt(x) <- -k * x\nt <- 2", "t <-")
  expect_identical(err$line, 2L)
  expect_identical(err$column, 1L)
  expect_compile_error("d/dt(x) <- -k * x\ny <- 2 * linCmt()", "linCmt()", "name <- linCmt()")
  expect_compile_error("cl <- 1\ny <- linCmt()", "linCmt()", "v (or vc, v1)")
  expect_compile_error("v <- 1\ny <- linCmt()", "linCmt()", "cl or k")
  expect_compile_error("cl <- 1\nv <- 2\na <- linCmt()\nb <- linCmt()", "linCmt()")
  expect_compile_error("d/dt(x) <- -k * x\nf(y) <- 0.5", "f(y)", "d/dt(y)")
  expect_compile_error("d/dt(x) <- -k * x\nx <- 2", "x <-", "x(0)")
  expect_compile_error("y <- 1", "model", "d/dt(x)")
  expect_compile_error("d/dt(x) <- -k * x\nalag(x) <- x / 10", "alag/f/rate/dur/x(0)")
  expect_compile_error("d/dt(x) <- -k * x\nf(x) <- t", "alag/f/rate/dur/x(0)")
  expect_compile_error("d/dt(x) <- -k * x\n2 + 2", "+")
  expect_compile_error("d/dt(x) <- -k * x\ny <- 'a'", "\"a\"")
})

test_that("a local read before every path assigns it is refused", {
  err <- expect_compile_error(
    "if (k > 1) {\n  a <- 2\n}\nd/dt(x) <- -a * x", "a",
    "an assignment to `a` on every path before"
  )
  expect_identical(err$line, 4L)
  # Assigned on both branches, it is fine.
  expect_s3_class(pmx_ode_model(
    "if (k > 1) {\n  a <- 2\n} else {\n  a <- 3\n}\nd/dt(x) <- -a * x"
  ), "pmx_ode_model")
  expect_compile_error("d/dt(x) <- -a * x\na <- 1", "a")
})

test_that("the corpus's one model outside the subset is refused for ifelse()", {
  texts <- corpus_model_texts()
  refused <- character()
  for (name in names(texts)) {
    e <- tryCatch({
      pmx_ode_model(texts[[name]])
      NULL
    }, pmx_ode_error_compile = function(e) e)
    if (!is.null(e)) refused[[name]] <- e$construct
  }
  expect_identical(refused, c(lib_tgi_no_sat_expo = "ifelse"))
})
