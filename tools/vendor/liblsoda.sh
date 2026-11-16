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
# 7. No abort() in solsy().
perl -0pi -e 's/\t\tabort\(\);/\t\treturn 0;/' "$dest/solsy.c"
# 6. (void) for empty parameter lists.
perl -pi -e 's/lsoda_create_ctx\(\)/lsoda_create_ctx(void)/; s/lsoda_create_opt\(\)/lsoda_create_opt(void)/' "$dest/lsoda.h" "$dest/lsoda.c"
perl -pi -e 's/void printcm12\(\)/void printcm12(void)/' "$dest/cfode_static.c"
cat >> "$dest/common.h" <<'EOT'

/* zupmxode: no console output from C (tools/vendor/liblsoda.sh). */
#define LSODA_LOG(...) ((void) 0)
EOT
