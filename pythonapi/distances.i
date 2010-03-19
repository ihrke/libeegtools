%{
#include "definitions.h"
#include "distances.h"
%}
typedef double  (*VectorDistanceFunction)       (double*,double*,int,OptArgList*);
double** matrix_init(int N, int M){
    int i,j;
    double **d;
    d = (double**) malloc( N*sizeof(double*) );
    for( i=0; i<N; i++){
    d[i] = (double*) malloc( M*sizeof(double) );
    for( j=0; j<M; j++ ){
    d[i][j]=0.0;
    }
    }
    
    return d;
    }

/* convert numpy 2D matrix to double** */
%typemap(in) (const double **X, int  n, int p){
  int is_new_object;
  PyArrayObject* array=NULL;
  array = obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE,
                                                   &is_new_object);
  if (!array || !require_dimensions(array, 2) ) SWIG_fail;
  double *d = (double*) array_data(array);
  int i,n,p;
  n = array_size( array, 0);
  p = array_size( array, 1);

  $1 = (double**)malloc( n*sizeof(double*));
  for( i=0; i<n; i++ ){
	 $1[i] = (double*)malloc( p*sizeof(double));
	 memcpy( $1[i], array+i*p, p*sizeof(double) );
  }
  $2=n;
  $3=p;
 }
%typemap(freearg) (const double **X, int  n, int p){
  
 }

%typemap(argout) (int n, int p, double **D){
  $result=doubleptrptr_to_pyarray( $3, $1, $2 );
 }

double** vectordist_distmatrix( VectorDistanceFunction f, const double **X, 
										  int n, int p, double **D, 
										  ProgressBarFunction progress=NULL, 
										  OptArgList* optargs=NULL );
%pythoncode 
%{ 
def distmatrix( X, distfct ):
   (d1,d2)=X.shape;
   D = buildMatrix( d1, d1 );
   return _pyeegtools.vectordist_distmatrix( distfct, X, D);
%}


/* this typemap converts a numpy array to a double pointer without any size checks
	-> use carefully! */
%typemap(in) (double *x1 ){  
  int is_new_object;
  PyArrayObject* array=NULL;
  array = obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE,
                                                   &is_new_object);
  if (!array || !require_dimensions(array, 1) ) SWIG_fail;
  $1 = (double*) array_data(array);
}
  //%apply (double* IN_ARRAY1[ANY]) { (double *) };
%apply (double* IN_ARRAY1, int DIM1){ (double *x2, int p)};

%callback("%s_cb");
/* these functions are also defined as callback functions */
double   vectordist_euclidean           ( double *x1, double *x2, int p, OptArgList *optargs=NULL ); 
double   vectordist_euclidean_normalized( double *x1, double *x2, int p, OptArgList *optargs=NULL );
double   vectordist_dtw                 ( double *x1, double *x2, int p, OptArgList *optargs=NULL );
double   vectordist_regularized_dtw     ( double *x1, double *x2, int p, OptArgList *optargs=NULL );
%nocallback;
