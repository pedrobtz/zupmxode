# Classed conditions (design §8). Every error the package raises on
# purpose has class `pmx_ode_error` and a subclass naming what failed, so
# callers and tests catch it by class, never by message.

pmx_abort <- function(subclass, message, ..., call = NULL) {
  cnd <- structure(
    class = c(subclass, "pmx_ode_error", "error", "condition"),
    list(message = message, call = call, ...)
  )
  stop(cnd)
}

pmx_warn <- function(subclass, message, ..., call = NULL) {
  cnd <- structure(
    class = c(subclass, "pmx_ode_warning", "warning", "condition"),
    list(message = message, call = call, ...)
  )
  warning(cnd)
}

# A verifier verdict (src/zpx_model.c's status list) as an error naming the
# block and instruction, 1-based as R counts.
abort_ir <- function(v, call = NULL) {
  where <- if (v$block >= 0) {
    paste0(" in block `", ir_blocks[v$block + 1L], "`",
           if (v$insn >= 0) paste0(", instruction ", v$insn + 1L) else "")
  } else {
    ""
  }
  pmx_abort("pmx_ode_error_ir", paste0("Invalid model IR", where, ": ",
                                       v$message, "."),
            block = if (v$block >= 0) ir_blocks[v$block + 1L] else NA_character_,
            insn = if (v$insn >= 0) v$insn + 1L else NA_integer_)
}
