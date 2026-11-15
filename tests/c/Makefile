# The engine's C tests, built without R: no R include path is given, so a
# file under src/core/ that includes R.h fails to compile here (roadmap
# Stage 0). Run from the package root: make -C tests/c
CC ?= cc
CFLAGS ?= -O2
STRICT = -std=c99 -Wall -Wextra -Wpedantic -Werror
CORE = ../../src/core
SRC = $(wildcard $(CORE)/*.c)

.PHONY: check clean
check: test_core
	./test_core

test_core: test_core.c $(SRC) $(wildcard $(CORE)/*.h)
	$(CC) $(STRICT) $(CFLAGS) -I$(CORE) -o $@ test_core.c $(SRC) -lm

clean:
	rm -f test_core test_core.exe
