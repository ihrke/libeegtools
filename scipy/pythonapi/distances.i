
%typemap(in) (double *x1, double *x2, int n){
  $1 = python_list_to_doubleptr($input);
  $2 = PyList_Size($input);
 }

/* Some callback functions */
%callback("%s_cb");
double   vectordist_euclidean           ( double *x1, double *x2, int p, void *userdata );
double   vectordist_euclidean_normalized( double *x1, double *x2, int p, void *userdata );
double   vectordist_dtw                 ( double *x1, double *x2, int p, void *userdata );
double   vectordist_regularized_dtw     ( double *x1, double *x2, int p, void *userdata );
%nocallback;
