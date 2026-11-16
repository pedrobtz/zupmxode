int stoda(struct lsoda_context_t * ctx, double *y, int jstart);
int      correction(struct lsoda_context_t * ctx, double *y, double pnorm, double *del, double *delp, double told, int *m);
int prja(struct lsoda_context_t * ctx, double *y);

int solsy(struct lsoda_context_t * ctx, double *y);
int      orderswitch(struct lsoda_context_t * ctx, double rhup, double dsm, double *rh, int kflag, int maxord);
int      intdy(struct lsoda_context_t * ctx, double t, int k, double *dky);
int      corfailure(struct lsoda_context_t * ctx, double told);
void     methodswitch(struct lsoda_context_t * ctx, double dsm, double pnorm, double *rh);
void     cfode(struct lsoda_context_t * ctx, int meth);
void     cfode_static(struct lsoda_context_t * ctx, int meth);
void     scaleh(struct lsoda_context_t * ctx, double rh);

#if defined(__GNUC__)
char * _strdup_printf(const char * fmt, ...) __attribute__((format(printf, 1, 2)));
#else
char * _strdup_printf(const char * fmt, ...);
#endif
#ifdef CFODE_STATIC
	#define cfode cfode_static
#endif
#define ERROR(...) (ctx->error? free(ctx->error):(void)1, ctx->error=_strdup_printf(__VA_ARGS__))
