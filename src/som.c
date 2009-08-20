/***************************************************************************
 *   Copyright (C) 2008/2009 by Matthias Ihrke   *
 *   mihrke@uni-goettingen.de   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "som.h"
#include "distances.h"
#include "mathadd.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


/** allocate the Som-struct 
 */
Som* som_init( int dimension, int n, int nruns, SOMConnectivityType connectivity ){
  Som *s;
  int i;
  s=(Som*)malloc(sizeof(Som));
  s->dimension=dimension;
  s->n=n;
  s->nruns = nruns;
  dprintf("dimension=%i, n=%i\n", s->dimension, s->n );
  s->m = (double**) malloc( n*sizeof(double*) );
  for( i=0; i<n; i++ ){
	 s->m[i] = (double*)malloc( dimension*sizeof(double));
  }
  s->distancefct = vectordist_euclidean;
  s->initial_runs = 0.1*nruns;
  s->neighbourhoodfct = som_neighbourhood_gaussian;
  s->connectivity = connectivity;
  if( s->connectivity==CUSTOM ){
	 s->custom_connectivity = (double**) malloc( n*sizeof(double*) );
	 for( i=0; i<n; i++ )
		s->custom_connectivity[i] = (double*) malloc( n*sizeof(double) );
  } 
  s->distancefct_parameters=NULL;  

  gsl_rng_env_setup();
  s->random_number_type = gsl_rng_default;
  s->rng = gsl_rng_alloc (s->random_number_type);

  return s;
}

void   som_free( Som *s ){
  int i;
  for( i=0; i<s->n; i++ ){
	 free( s->m[i] );
  }
  free( s->m );
  if(  s->connectivity==CUSTOM ){
	 for( i=0; i<s->n; i++ )
		free( s->custom_connectivity[i] );
	 free( s->custom_connectivity ) ;
  }
  free( s );
}

double som_neighbourhood_gaussian( int x, int bmu, struct som_struct *s, int t ){
  double sigma, alpha;	  
  double h;
  /* sigma shrinks with runs, fast during initial phase, then slower */
  sigma = exp(-(double)t/(double)s->initial_runs)*(s->n/2) + 1;
  
  if(t<=s->initial_runs)
	 alpha = 0.9*(1-(double)t/(double)(s->initial_runs));
  else
	 alpha = 0.02*(1-(double)t/(double)s->nruns);
  
  if( s->connectivity==ONED_LINEAR ){
	 h = alpha * exp(-pow( ABS(x-bmu), 2.0)/(double)(2.0*pow(sigma, 2.0)));
  } else if( s->connectivity!=ONED_LINEAR ) {
	 errprintf("Other topologies than ONED_LINEAR are not implemented yet!\n");
	 return -1;
  }
  return h;
}

/** apply SOM-training rule
	 \f[
	 m_i(t+1) = m_i(t) + h_{ci}(t)[ x(t) - m_i(t) ]
	 \f]
	 
 */
void som_train_from_data( Som *s, double **X, int dim, int nsamples ){
  int i, j, t;
  int bmu; /* best matching unit */
  double h, tmp; 
  double bmu_score; /* best matching unit score */
  double *input;

  if( dim != s->dimension ){
	 errprintf( "Data dimension and SOM-dimension do not match... exit\n");
	 return;
  }

  if( s->progress ){
	 (*(s->progress))( PROGRESSBAR_INIT, s->nruns );
  }

  for( t=1; t<=s->nruns; t++ ){
	 if( s->progress ){
		s->progress( PROGRESSBAR_CONTINUE_LONG, t-1 );
	 }

	 /* choose a sample */
	 input = X[(t-1) % nsamples];

	 /* find BMU */
	 for( i=0; i<s->n; i++ ){	 
		if( (tmp=s->distancefct( input, s->m[i], dim, NULL ))<bmu_score ){
		  bmu_score = tmp;
		  bmu = i;
		}
	 }
	 
	 /* adapt weights (learning) */
	 for( i=0; i<s->n; i++ ){
		if( s->progress ){
		  s->progress( PROGRESSBAR_CONTINUE_SHORT, i );
		}

		h = s->neighbourhoodfct( i, bmu, s, t );
		for( j=0; j<dim; j++ ){
		  s->m[i][j] += h*(input[j]-s->m[i][j]);
		}
	 }
  }

}


/** initialize the codebook vectors with random values between 
	 min and max (uniformely drawn);
*/
void som_initialize_random( Som *s, double min, double max ){
  int i,j;

  for( i=0; i<s->n; i++ ){
	 for( j=0; j<s->dimension; j++ ){
		s->m[i][j] = gsl_ran_flat( s->rng, min, max );
	 }
  }
}


/** initialize the codebook vectors with randomly drawn
	 samples from the data X; 

	 \param s the model
	 \param X data of dimension [dim, nsamples]
*/
void som_initialize_random_samples( Som *s, double **X, int dim, int nsamples ){
  int i;
  int ransamp;
  
  if( dim != s->dimension ){
	 errprintf( "Data dimension and SOM-dimension do not match... exit\n");
	 return;
  }

  for( i=0; i<s->n; i++ ){
	 ransamp = gsl_rng_uniform_int( s->rng, nsamples );
	 memcpy( s->m[i], X[ransamp], dim*sizeof(double) );
  }
}

void som_print( FILE *f, Som *s ){
  fprintf( f, 
			  "Som-struct '%p'\n"
			  "  dimension = %i\n"
			  "  n         = %i\n"
			  "  nruns     = %i\n"
			  "  initial   = %i\n", s, s->dimension, s->n, s->nruns, 
			  s->initial_runs );
}