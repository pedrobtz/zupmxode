# The engine's C tests, built without R: no R include path is given, so a
# file under src/core/ that includes R.h fails to compile here (roadmap
# Stage 0). Run from the package root: make -C tests/c
CC ?= cc
CFLAGS ?= -O2
STRICT = -std=c99 -Wall -Wextra -Wpedantic -Werror
CORE = ../../src/core
SRC = $(wildcard $(CORE)/*.c)
LSODA = ../../src/vendor/liblsoda
LSODA_SRC = $(wildcard $(LSODA)/*.c)

.PHONY: check clean
check: test_core test_lsoda
	./test_core
	./test_lsoda

# The vendored library is not held to -Werror -Wpedantic: it uses GNU's
# , ## __VA_ARGS__, and is upstream's code but for tools/vendor's patches.
test_lsoda: test_lsoda.c $(LSODA_SRC) $(wildcard $(LSODA)/*.h)
	$(CC) -std=gnu99 $(CFLAGS) -I$(LSODA) -o $@ test_lsoda.c $(LSODA_SRC) -lm

test_core: test_core.c $(SRC) $(wildcard $(CORE)/*.h)
	$(CC) $(STRICT) $(CFLAGS) -I$(CORE) -o $@ test_core.c $(SRC) -lm

clean:
	rm -f test_core test_core.exe test_lsoda test_lsoda.exe
