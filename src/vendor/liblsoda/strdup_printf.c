#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
char * _strdup_printf(const char * fmt, ...) {
	char one[1];
	if (fmt == NULL) return NULL;
	va_list va;
	va_start(va, fmt);
	int s = vsnprintf(one, sizeof one, fmt, va);
	va_end(va);
	char * rt = malloc((size_t) s + 1);
	if (rt == NULL) return NULL;
	va_start(va, fmt);
	vsnprintf(rt, (size_t) s + 1, fmt, va);
	va_end(va);
	return rt;
}
