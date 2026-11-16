#!/bin/sh
# Vendors liblsoda (design D5) into src/vendor/liblsoda/ from the commit in
# tools/vendor/manifest.tsv, then applies zupmxode's patches, which are all
# here and nowhere else:
#
#   1. Every printf() and fprintf(stderr, ...) becomes LSODA_LOG(...), a
#      macro that does nothing (common.h): CRAN forbids C code that writes
#      to the console, and the messages are in ctx->error anyway.
#   2. Only the sources upstream's src/Makefile builds are taken (not
#      printcf.c, debugging output, nor ewset.c, an #error stub).
#   4. The six BLAS routines' K&R definitions become the prototypes blas.h
#      declares: C23, which R-devel's compilers default to, has no K&R.
#   5. lsoda(): ihit starts at 0, and softfailure()'s loop variable is i_,
#      not i, which shadowed the caller's i that its message prints (an
#      uninitialised read, -Wuninitialized).
#   6. Functions declared with empty parentheses take (void).
#   7. solsy() returns 0 instead of calling abort() for a banded Jacobian,
#      which liblsoda never builds (miter is always 2): CRAN forbids C that
#      can end the R process.
#   8. lsoda(): a stray comma in the message for a repeatedly failed error
#      test made the rest of the string literal into arguments, so %g read
#      a char * and the appended %s a double: a crash on that path (found by
#      ASan in zupmxode's Stage 6).
#   9. ERROR(), hardfailure() and softfailure() take their format in
#      __VA_ARGS__, portable C99, not GNU's , ## __VA_ARGS__; the message's
#      " @(file:line)" suffix goes with it.
#  11. The BLAS routines are renamed lsoda_daxpy, lsoda_ddot, lsoda_dgefa,
#      lsoda_dgesl, lsoda_dscal, lsoda_idamax: under the plain names the
#      dynamic linker bound liblsoda's calls to the BLAS R links (OpenBLAS's
#      Fortran idamax_, which takes pointers), reading an int as an address
#      (found by gcc ASan, Stage 6).
#  10. _strdup_printf() is declared with printf's format attribute, so the
#      compiler checks every message against its arguments, and measures
#      the length into a one-byte buffer rather than NULL, returning NULL
#      for a NULL format (gcc's fortified headers, under UBSan, warn of a
#      null format on that path).
#  12. lsoda_prepare() returns 0 when its calloc() fails, and alloc_mem()
#      when its malloc() does, before either writes through the pointer:
#      upstream dereferenced both unchecked (found by zupmxode's
#      alloc-failure sweep, Stage 12).
#   3. _strdup_printf() allocates one byte more, for the string's NUL, and
#      returns NULL when malloc() does: upstream truncated every message by
#      a character and wrote through an unchecked pointer.
#
#   sh tools/vendor/liblsoda.sh [destination]    # default src/vendor/liblsoda
set -eu
dest=${1:-src/vendor/liblsoda}
commit=$(awk -F '\t' '$1 == "liblsoda" { print $3 }' tools/vendor/manifest.tsv)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
curl -fsSL "https://codeload.github.com/sdwfrost/liblsoda/tar.gz/$commit" -o "$tmp/l.tgz"
tar xzf "$tmp/l.tgz" -C "$tmp"
src="$tmp/liblsoda-$commit/src"
rm -rf "$dest"
mkdir -p "$dest"
srcs=$(sed -n '/^SRCS/,/^$/p' "$src/Makefile" | tr -d '\\' | tr ' ' '\n' | grep '\.c$')
for f in $srcs; do cp "$src/$f" "$dest/"; done
cp "$src"/*.h "$src"/cfode_static.inc "$dest/"
cp "$tmp/liblsoda-$commit/LICENSE" "$dest/LICENSE"
for f in "$dest"/*.c; do
  # Portable sed: no \b in BSD sed.
  sed -e 's/fprintf(stderr, */LSODA_LOG(/g' \
      -e 's/\([^A-Za-z0-9_]\)printf(/\1LSODA_LOG(/g' -e 's/^printf(/LSODA_LOG(/' "$f" > "$f.new"
  mv "$f.new" "$f"
done
sed -e 's/char \* rt = malloc(s);/char * rt = malloc((size_t) s + 1);\
	if (rt == NULL) return NULL;/' -e 's/vsnprintf(rt, s, fmt, va);/vsnprintf(rt, (size_t) s + 1, fmt, va);/' \
  "$dest/strdup_printf.c" > "$dest/strdup_printf.c.new"
mv "$dest/strdup_printf.c.new" "$dest/strdup_printf.c"
# 4. K&R to prototypes (perl -0: each definition spans three lines).
perl -0pi -e 's/daxpy\(n, da, dx, incx, dy, incy\)\n\tdouble +da, \*dx, \*dy;\n\tint +n, incx, incy;\n/daxpy(int n, double da, double dx[], int incx, double dy[], int incy)\n/' "$dest/daxpy.c"
perl -0pi -e 's/ddot\(n, dx, incx, dy, incy\)\n\tdouble +\*dx, \*dy;\n\tint +n, incx, incy;\n/ddot(int n, double dx[], int incx, double dy[], int incy)\n/' "$dest/ddot.c"
perl -0pi -e 's/dgefa\(a, n, ipvt, info\)\n\tdouble +\*\*a;\n\tint +n, \*ipvt, \*info;\n/dgefa(double **a, int n, int *ipvt, int *info)\n/' "$dest/dgefa.c"
perl -0pi -e 's/dgesl\(a, n, ipvt, b, job\)\n\tdouble +\*\*a, \*b;\n\tint +n, \*ipvt, job;\n/dgesl(double **a, int n, int *ipvt, double b[], int job)\n/' "$dest/dgesl.c"
perl -0pi -e 's/dscal\(n, da, dx, incx\)\n\tdouble +da, \*dx;\n\tint +n, incx;\n/dscal(int n, double da, double dx[], int incx)\n/' "$dest/dscal.c"
perl -0pi -e 's/idamax\(n, dx, incx\)\n\tdouble +\*dx;\n\tint +n, incx;\n/idamax(int n, double dx[], int incx)\n/' "$dest/idamax.c"
# 5. Uninitialised reads in lsoda.c.
perl -0pi -e 's/int             i=0, ihit;/int             i=0, ihit = 0;/' "$dest/lsoda.c"
perl -0pi -e 's/(#define softfailure\(code, fmt,\.\.\.\) \\\n\{ \\\n)\tint             i; \\\n/$1\tint             i_; \\\n/; s/\tfor \(i = 1; i <= neq; i\+\+\) \\\n\t  y\[i\] = _C\(yh\)\[1\]\[i\];/\tfor (i_ = 1; i_ <= neq; i_++) \\\n\t  y[i_] = _C(yh)[1][i_];/' "$dest/lsoda.c"
# 8. The stray comma.
perl -0pi -e 's/(softfailure\(-4, "lsoda -- at t = %g and step size _C\(h\) = %g, the\\n"),\n/$1\n/' "$dest/lsoda.c"
# 9. C99 variadic macros.
perl -0pi -e 's/#define ERROR\(fmt, \.\.\.\) \(ctx->error\? free\(ctx->error\):\(void\)1, ctx->error=_strdup_printf\("EE:" fmt " @\(%s:%d\)", ## __VA_ARGS__, __FILE__, __LINE__\)\)/#define ERROR(...) (ctx->error? free(ctx->error):(void)1, ctx->error=_strdup_printf(__VA_ARGS__))/' "$dest/lsoda_internal.h"
perl -0pi -e 's/#define hardfailure\(fmt,\.\.\.\) \\\n\{ \\\n\tERROR\(fmt, ## __VA_ARGS__\);/#define hardfailure(...) \\\n{ \\\n\tERROR(__VA_ARGS__);/; s/#define softfailure\(code, fmt,\.\.\.\) \\\n(\{ \\\n\tint             i_; \\\n\tint neq = ctx->neq; \\\n \\\n)\tERROR\(fmt, ## __VA_ARGS__\);/#define softfailure(code, ...) \\\n$1\tERROR(__VA_ARGS__);/' "$dest/lsoda.c"
# 10. The format attribute, and no NULL buffer.
perl -0pi -e 's/char \* _strdup_printf\(char \* fmt, \.\.\.\);/#if defined(__GNUC__)\nchar * _strdup_printf(const char * fmt, ...) __attribute__((format(printf, 1, 2)));\n#else\nchar * _strdup_printf(const char * fmt, ...);\n#endif/' "$dest/lsoda_internal.h"
perl -0pi -e 's/char \* _strdup_printf\(char \* fmt, \.\.\.\) \{/char * _strdup_printf(const char * fmt, ...) {\n\tchar one[1];\n\tif (fmt == NULL) return NULL;/; s/int s = vsnprintf\(NULL, 0, fmt, va\);/int s = vsnprintf(one, sizeof one, fmt, va);/' "$dest/strdup_printf.c"
# 11. Rename the BLAS routines, at every use and definition.
for f in "$dest"/*.c "$dest"/*.h; do
  perl -pi -e 's/\b(daxpy|ddot|dgefa|dgesl|dscal|idamax)\b/lsoda_$1/g' "$f"
done
# 12. Check the two allocations before they are used.
perl -0pi -e 's/(\tctx->common = calloc\(1, sizeof\(struct lsoda_common_t\)\);\n)/$1\tif (ctx->common == NULL) return 0;\n/; s/(\t_C\(memory\) = malloc\(offset\);\n)/$1\tif (_C(memory) == NULL) return 0;\n/' "$dest/lsoda.c"
grep -q 'if (ctx->common == NULL) return 0;' "$dest/lsoda.c"
grep -q 'if (_C(memory) == NULL) return 0;' "$dest/lsoda.c"
# 7. No abort() in solsy().
perl -0pi -e 's/\t\tabort\(\);/\t\treturn 0;/' "$dest/solsy.c"
# 6. (void) for empty parameter lists.
perl -pi -e 's/lsoda_create_ctx\(\)/lsoda_create_ctx(void)/; s/lsoda_create_opt\(\)/lsoda_create_opt(void)/' "$dest/lsoda.h" "$dest/lsoda.c"
perl -pi -e 's/void printcm12\(\)/void printcm12(void)/' "$dest/cfode_static.c"
cat >> "$dest/common.h" <<'EOT'

/* zupmxode: no console output from C (tools/vendor/liblsoda.sh). */
#define LSODA_LOG(...) ((void) 0)
EOT
