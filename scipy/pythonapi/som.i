  
  /**\ingroup som_help 
	\{ */
  Som* som_init( int dimension, int n, int nruns, SOMConnectivityType connectivity );
  void som_free( Som *s );
  void som_print( FILE *f, Som *s );
  /** \} */

  /**\ingroup som_initialize
	\{ */ 
  void som_initialize_random( Som *s, double min, double max );
  void som_initialize_random_samples( Som *s, double **X, int dim, int nsamples );
  double** som_generate_connectivity_matrix( SOMConnectivityType type, double **m, int n );
  /** \} */

  /**\ingroup som_neighbourhood
	\{ */
  double som_neighbourhood_gaussian( int x, int bmu, struct som_struct *s, int t);
  /** \} */

  double som_time_decay_linear( int t, int nruns, int initial_runs );

  /**\ingroup som_train
	\{ */
  void som_train_from_data( Som *s, double **X, int dim, int nsamples );
  /** \} */
  
