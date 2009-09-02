%{
#include <libeegtools/mathadd.h>
%}
/*%include "numpy.i" 

%init %{ 
import_array(); 
%} 
*/


#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <libeegtools/helper.h>
#include <libeegtools/definitions.h>


double  mad(const double *data, int n); 
int     sgn(int x);
double  maxel(double *v, int n);
int     maxeli(int *v, int n);
int     closest_index(const double *v, int n, double c);

double* sampled_line(double *ntimes, int n, double start, double end);
double* lininterp(const double *x1, const double *y1, int n1, 
						const double *x2,       double *y2, int n2);
int*    linspace(int first, int last);

int     bresenham_howmany_points( int xstart,int ystart,int xend,int yend );
int*    bresenham(int xstart,int ystart,int xend,int yend, int *points);

double  gaussian( double x, double sigma, double mu );

double** disttransform_deadreckoning(int **I, int X, int Y, double **d);

int     next_pow2( int n );
int     iremainder( double x, double y);
Complex* expand_polynomial_from_roots( const Complex *roots, int n, Complex *coeffs );


void fft(double *data,  unsigned long nn, int isign);
  

double  drawsample_nearest_neighbour( const double *v, int n, double x );
double* resample_linear( const double *s, int n, int newn, double *news );
double* resample_nearest_neighbour( const double *s, int n, int newn, double *news );
double* resample_gsl( const double *s, int n, int newn, double *news, const gsl_interp_type *method );


double rmse(const double *r, const double *d, int n);
double snr (const double *r, const double *d, int n);


double  vector_min( double *v, int n, int *idx );
double* vector_init( double *v, int n,  double val );
void    vector_minus_scalar( double *v, int n, double val );
double* vector_complex_to_real( const Complex *vc, double *vr, int n );
double  vector_euclidean_distance( const double *v1, const double *v2, int n );
void    vector_shuffle_int( int *permut, int n );


double** matrix_delrow(double **m, int N, int n, int row);
double** matrix_delcol(double **m, int N, int n, int col);
double   matrix_min(const double **m, int N, int n, int *i1, int *i2);
double   matrix_max(const double **m, int N, int n, int *i1, int *i2);
void     matrix_print(double **m, int N, int n);
double** matrix_init(int N, int M);
int**    matrix_init_int(int N, int M);
void     matrix_divide_scalar(double **m, int N, int n, double s);
void     matrix_add_scalar(double **m, int N, int n, double s);
void     matrix_mul_scalar(double **m, int N, int n, double s);
void     matrix_normalize_by_max( double **m, int M, int N );
void     matrix_add_matrix(double **m1, const double **m2, int N, int n);
void     matrix_dottimes_matrix( double **m1, const double **m2, int N, int M );
void     matrix_copy( const double **src, double **dest, int N, int M );
void     scalar_minus_matrix( double scalar, double **m, int N, int M );
double** matrix_rand( double **m, int N, int M, double lower, double upper );
void     matrix_free(double **m, int N);

