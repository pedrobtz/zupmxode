# The IR (design §4.1, D13), on the R side: an integer vector of header and
# code plus a constant pool, laid out as src/core/ir.h describes. The C
# verifier is the authority on what is valid; nothing here is trusted by C.

ir_magic <- 0x5a50584fL
ir_version <- 1L
ir_blocks <- c("pre", "ivl", "rhs", "post", "err")
ir_width <- 5L

# Opcodes and their operand kinds, in the order the IR stores them: r a
# register, k a constant, i a slot index, f a math function, j a label.
ir_ops <- list(
  const = list(code = 1L, args = c("r", "k")),
  loadt = list(code = 2L, args = "r"),
  loady = list(code = 3L, args = c("r", "i")),
  loadp = list(code = 4L, args = c("r", "i")),
  loadw = list(code = 5L, args = c("r", "i")),
  storeo = list(code = 6L, args = c("i", "r")),
  storew = list(code = 7L, args = c("i", "r")),
  mov = list(code = 8L, args = c("r", "r")),
  add = list(code = 10L, args = c("r", "r", "r")),
  sub = list(code = 11L, args = c("r", "r", "r")),
  mul = list(code = 12L, args = c("r", "r", "r")),
  div = list(code = 13L, args = c("r", "r", "r")),
  neg = list(code = 14L, args = c("r", "r")),
  call1 = list(code = 20L, args = c("r", "f", "r")),
  call2 = list(code = 21L, args = c("r", "f", "r", "r")),
  lt = list(code = 30L, args = c("r", "r", "r")),
  le = list(code = 31L, args = c("r", "r", "r")),
  gt = list(code = 32L, args = c("r", "r", "r")),
  ge = list(code = 33L, args = c("r", "r", "r")),
  eq = list(code = 34L, args = c("r", "r", "r")),
  ne = list(code = 35L, args = c("r", "r", "r")),
  and = list(code = 36L, args = c("r", "r", "r")),
  or = list(code = 37L, args = c("r", "r", "r")),
  not = list(code = 38L, args = c("r", "r")),
  jmp = list(code = 40L, args = "j"),
  jmpz = list(code = 41L, args = c("r", "j"))
)

# The math table, by the names the model language uses (design §4.1); the
# numbers are src/core/ir.h's.
ir_math <- c(
  exp = 0L, log = 1L, sqrt = 2L, abs = 3L, floor = 4L, ceiling = 5L,
  sin = 6L, cos = 7L, tan = 8L, gamma = 9L, lgamma = 10L, logit = 11L,
  expit = 12L, probit = 13L, probitInv = 14L,
  pow = 32L, min = 33L, max = 34L
)

# Assembles IR from text, one instruction per line: `op operand ...`, with
# registers as r0, r1, ..., constants written as numbers (pooled), slots as
# integers, math functions by name, and labels as `name:` lines that jumps
# name. `blocks` is a named list of texts, by block name; `nout` the output
# length of each block given. For tests and, from Stage 3, the compiler.
ir_asm <- function(blocks = list(), ny = 0L, npar = 0L, nwork = 0L,
                   nout = integer()) {
  stopifnot(all(names(blocks) %in% ir_blocks),
            all(names(nout) %in% ir_blocks))
  consts <- double()
  pool <- function(x) {
    i <- match(x, consts)
    # match() finds NaN and -0 by value; the pool keeps them distinct.
    if (is.na(i) || !identical(consts[i], x, num.eq = FALSE)) {
      consts[[length(consts) + 1L]] <<- x
      i <- length(consts)
    }
    i - 1L
  }
  desc <- integer()
  code <- integer()
  for (b in ir_blocks) {
    lines <- trimws(strsplit(if (is.null(blocks[[b]])) "" else blocks[[b]],
                             "\n", fixed = TRUE)[[1]])
    lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
    is_label <- grepl(":$", lines)
    labels <- stats::setNames(cumsum(!is_label)[is_label], sub(":$", "", lines[is_label]))
    insns <- lines[!is_label]
    nreg <- 0L
    words <- strsplit(insns, "[[:space:]]+")
    rows <- vapply(words, function(w) {
      op <- ir_ops[[w[1]]]
      if (is.null(op)) stop("unknown IR opcode `", w[1], "`", call. = FALSE)
      args <- w[-1]
      if (length(args) != length(op$args)) {
        stop("`", w[1], "` takes ", length(op$args), " operands", call. = FALSE)
      }
      vals <- vapply(seq_along(args), function(k) {
        a <- args[k]
        switch(op$args[k],
          r = {
            r <- as.integer(sub("^r", "", a))
            nreg <<- max(nreg, r + 1L)
            r
          },
          k = pool(ir_number(a)),
          i = as.integer(a),
          f = unname(ir_math[a]),
          j = unname(labels[a])
        )
      }, 1L)
      if (anyNA(vals)) stop("bad operand in `", paste(w, collapse = " "), "`",
                            call. = FALSE)
      c(op$code, vals, integer(ir_width - 1L - length(vals)))
    }, integer(ir_width))
    out <- if (is.null(nout[b]) || is.na(nout[b])) 0L else as.integer(nout[[b]])
    desc <- c(desc, nreg, out, length(insns))
    code <- c(code, as.vector(rows))
  }
  new_ir(c(ir_magic, ir_version, as.integer(c(ny, npar, nwork)),
           length(consts), desc, code), consts)
}

# A constant as the assembler reads it: "b" and the 16 hexadecimal digits of
# its IEEE bits, as the compiler writes every constant, or a decimal
# number, as tests write them by hand. The bits are what the compiler
# uses: R's decimal reader is exact only where long double has more bits
# than double, and is not under valgrind or on macOS arm64.
ir_number <- function(a) {
  if (grepl("^b[0-9a-f]{16}$", a)) {
    bytes <- as.raw(strtoi(substring(a, seq(2, 16, 2), seq(3, 17, 2)), 16L))
    return(readBin(bytes, "double", size = 8, endian = "big"))
  }
  as.numeric(a)
}

# The constant's assembler form: "b" and its IEEE bits.
ir_bits <- function(x) {
  paste0("b", paste(format(writeBin(as.double(x), raw(), endian = "big")), collapse = ""))
}

new_ir <- function(code, consts) {
  structure(list(code = as.integer(code), consts = as.double(consts)),
            class = "zupmxode_ir")
}

# The verifier's verdict: list(status, message, block, insn), status "ok"
# or "ir"; block and insn are 0-based, -1 when not at one.
ir_verify <- function(ir) {
  .Call(zpx_ir_verify_r, ir$code, ir$consts)
}

ir_check <- function(ir, call = parent.frame()) {
  v <- ir_verify(ir)
  if (v$status != "ok") abort_ir(v, call)
  invisible(ir)
}

ir_dims <- function(ir) {
  h <- ir$code
  list(ny = h[3], npar = h[4], nwork = h[5],
       nout = stats::setNames(h[7 + 3L * (seq_along(ir_blocks) - 1L) + 1L], ir_blocks))
}

# Block b by the reference evaluator (src/core/ir.c): list(out, work).
ir_eval <- function(ir, block, t = 0, y = double(), par = double(),
                    work = NULL) {
  d <- ir_dims(ir)
  if (is.null(work)) work <- double(d$nwork)
  .Call(zpx_ir_eval_r, ir$code, ir$consts, match(block, ir_blocks) - 1L,
        as.double(t), as.double(y), as.double(par), as.double(work))
}
