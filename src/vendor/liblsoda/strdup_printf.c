#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
char * _strdup_printf(char * fmt, ...) {
	va_list va;
	va_start(va, fmt);
	int s = vsnprintf(NULL, 0, fmt, va);
	va_end(va);
	char * rt = malloc((size_t) s + 1);
	if (rt == NULL) return NULL;
	va_start(va, fmt);
	vsnprintf(rt, (size_t) s + 1, fmt, va);
	va_end(va);
	return rt;
}
