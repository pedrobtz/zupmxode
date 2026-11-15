/* Builds src/core/ with no R headers on the include path and checks what
 * it reports. Exit status 0 on success. */
#include <stdio.h>
#include "core.h"

int main(void) {
  int failed = 0;
  if (zupmxode_core_abi() != ZUPMXODE_CORE_ABI) {
    fprintf(stderr, "core ABI %d, header says %d\n", zupmxode_core_abi(),
            ZUPMXODE_CORE_ABI);
    failed = 1;
  }
  if (!failed) printf("core: ok (ABI %d)\n", ZUPMXODE_CORE_ABI);
  return failed;
}
