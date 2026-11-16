
double lsoda_ddot(int n, double dx[], int incx, double dy[], int incy);
void lsoda_dgesl(double **a, int n, int * ipvt, double b[], int job);
void lsoda_dgefa(double ** a, int n, int * ipvt, int * info);
void lsoda_daxpy(int n, double da, double dx[], int incx, double dy[], int incy);
int lsoda_idamax(int n, double dx[], int incx);
void lsoda_dscal(int n, double da, double dx[], int incx);
double vmnorm(int n, double *v, double *w);
double fnorm(int n, double **a, double *w);
#if 0
static double vmnorm(int n, double *v, double *w)

/*
   This function routine computes the weighted max-norm
   of the vector of length n contained in the array v, with weights
   contained in the array w of length n.

   vmnorm = max( i = 1, ..., n ) fabs( v[i] ) * w[i].
*/

{
	int             i;
	double          vm;

	vm = 0.;
	for (i = 1; i <= n; i++)
		vm = fmax(vm, fabs(v[i]) * w[i]);
	return vm;

}
#endif
