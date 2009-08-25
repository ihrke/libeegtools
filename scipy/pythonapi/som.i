 
typedef double(*NeighbourhoodFunction)(int,int,struct som_struct*, int t);
typedef enum {
  ONED_LINEAR,
  TWOD_GRID,
  TWOD_HEXAGONAL,
  CUSTOM
} SOMConnectivityType;

typedef struct som_struct {
  double  **m; /**< code-book vectors */
  int     dimension; /**< dimension of code-book vector (corresponds to data-dimensionality) */
  int     n; /**< how many code-book vectors? */
  VectorDistanceFunction distancefct; /**< distance measure between data and codebook-vectors */
  void    *distancefct_parameters; /**< optional; pointer to parameters for some 
												  of the VectorDistanceFunction s */
  
  int     nruns; /**< number of runs to convergence */
  int     initial_runs; /**< number of runs with large flexibility (ordering phase) after
									which it is more restricted; e.g. 0.1*nruns */
  NeighbourhoodFunction neighbourhoodfct; /**< neighbourhood function  */
  SOMConnectivityType connectivity; /**< giving the dimension/structure of the SOM */
  double **custom_connectivity; /**< custom connectivity matrix that can implement
											  an arbitrary topology for the SOM-network */
  
  gsl_rng_type *random_number_type; /**< GSL-random number generator type */
  gsl_rng *rng; /**< the GSL random number generator */
  
  ProgressBarFunction progress; 
} Som;



Som* som_init( int dimension, int n, int nruns, SOMConnectivityType connectivity );
void som_free( Som *s );
void som_print( FILE *f, Som *s );

void som_initialize_random( Som *s, double min, double max );
void som_initialize_random_samples( Som *s, double **X, int dim, int nsamples );

double som_neighbourhood_gaussian( int x, int bmu, struct som_struct *s, int t);

void som_train_from_data( Som *s, double **X, int dim, int nsamples );
  
