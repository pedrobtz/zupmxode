/* inst/include/zupmxode.h compiled on its own, as C99 with every warning
 * an error and no R header (roadmap Stage 10): a consumer's estimator
 * can include it anywhere. Uses each type, so a typo in one fails here. */
#include "zupmxode.h"

#include <stdio.h>

int main(void) {
  zupmxode_eval_opts o = {ZUPMXODE_EVAL_OPTS_VERSION, 0, ZUPMXODE_ERR_COMBINED, 0.5, 0.1,
                          ZUPMXODE_METHOD_DEFAULT, 0, 0, 10};
  zupmxode_api api = {0};
  const zupmxode_model *m = NULL;
  const zupmxode_subject *s = NULL;
  zupmxode_arena *a = NULL;
  (void) m;
  (void) s;
  (void) a;
  int ok = api.subject_eval == NULL && o.version == 1 && ZUPMXODE_API_VERSION == 1 &&
           ZUPMXODE_E_FIRST_ERROR > ZUPMXODE_SS_UNCONVERGED && ZUPMXODE_OK == 0;
  printf("zupmxode.h standalone: %s\n", ok ? "ok" : "FAILED");
  return ok ? 0 : 1;
}
