/* The engine's core: C that builds without R's headers, so the C tests
 * under tests/c, the fuzz harnesses and mutation runs can compile it on
 * their own (roadmap Stage 0). Nothing under src/core/ may include R.h,
 * Rinternals.h or any header that does; src/zupmxode.c is the glue. */
#ifndef ZUPMXODE_CORE_H
#define ZUPMXODE_CORE_H

/* The core's ABI version, which the R glue reports and the C tests check.
 * Bumped when a struct or a function signature under src/core/ changes. */
#define ZUPMXODE_CORE_ABI 1

int zupmxode_core_abi(void);

#endif
