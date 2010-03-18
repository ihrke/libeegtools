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

/** generate a connectivity matrix of type type. 
	 The (i,j)th entry gives the distance between nodes i and j in
	 number of nodes.

	 type:
	 - ONED_LINEAR: \f$ m_{i,j} = | i-j | \f$
	 - TWOD_GRID: n should be a square number and \f$ m_{i,j} = \f$
	 - TWOD_HEXAGONAL: 

	 \param type 
	 \param m if NULL or ALLOC_IN_FCT, allocate own memory.
	 \param n number of nodes
 */
double** som_generate_connectivity_matrix( SOMConnectivityType type, double **m, int n ){
  int i,j;
  int sqrtn=0;
  int ix,iy;
  int jx,jy;	
  double diff;
  if( m==ALLOC_IN_FCT) {
	 m = dblpp_init( n,n );
  }
  switch( type ){
  case ONED_LINEAR:
	 for( i=0; i<n; i++ ){
		for( j=i; j<n; j++ ){
		  m[i][j] = ABS( i-j );
		  m[j][i] = m[i][j];
		}
	 }
	 break;
  case TWOD_GRID:
	 if( SQR( sqrtn=sqrt((double)n) ) != n ){
		warnprintf("for a 2D-grid, n should be the square of an integer\n");
	 }
	 for( i=0; i<n; i++ ){
		ix=i/sqrtn;
		iy=i-sqrtn*ix;
		for( j=i; j<n; j++ ){
		  jx = j/sqrtn;
		  jy = j-jx*sqrtn;
		  m[i][j] = sqrt( SQR(ix-jx)+SQR(iy-jy) );
		  m[j][i] = m[i][j];
		}
	 }	 
	 break; 
  case TWOD_HEXAGONAL:
	 for( i=0; i<n; i++ ){
		ix=i/sqrtn;
		iy=i-sqrtn*ix;
		for( j=i; j<n; j++ ){
		  jx = j/sqrtn;
		  jy = j-jx*sqrtn;
		  diff = ix-jx;
		  if( (iy-jy)%2 != 0 ){
			 if( iy%2 == 0 )
				diff -= 0.5;
			 else 
				diff += 0.5;
		  }
		  m[i][j] = SQR( diff );
		  diff = iy-jy;
		  m[i][j] += 0.75*SQR( diff );
		  m[i][j] = sqrt( m[i][j] );
		  m[j][i] = m[i][j];
		}
	 }	 
	 break; 
  default:
	 errprintf( "Do not know type '%i'\n", type );
  }
  return m;
}

/** allocate the Som-struct 
 */
Som* som_init( int dimension, int n, int nruns, SOMConnectivityType connectivity_type ){
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
  s->time_decay = som_time_decay_linear;
  s->initial_runs = 0.1*nruns;
  s->neighbourhoodfct = som_neighbourhood_gaussian;
  s->connectivity_type = connectivity_type;
  if( s->connectivity_type!=CUSTOM ){
	 s->connectivity = som_generate_connectivity_matrix( s->connectivity_type, 
																		  ALLOC_IN_FCT, s->n );
  } 
  s->distancefct_parameters=NULL;  

  gsl_rng_env_setup();
  s->random_number_type = (gsl_rng_type *)gsl_rng_default;
  s->rng = gsl_rng_alloc (s->random_number_type);

  return s;
}

void   som_free( Som *s ){
  int i;
  for( i=0; i<s->n; i++ ){
	 free( s->m[i] );
  }
  free( s->m );
  if(  s->connectivity_type==CUSTOM ){
	 for( i=0; i<s->n; i++ )
		free( s->connectivity[i] );
	 free( s->connectivity ) ;
  }
  free( s );
}


/** Linear time-decay function between 0 and 1.
	 \f[
	 f(t/\mbox{nruns}) = \begin{cases}
	 x & \mbox{if}\\
	 y & \mbox{otherwise}
	 \end{cases}
	 \f]

 */
double som_time_decay_linear( int t, int nruns, int initial_runs ){
  double alpha;

  if( t<=initial_runs ){
	 alpha = 0.9*(1-(double)t/(double)(initial_runs));
  } else {
	 alpha = 0.02*(1-(double)t/(double)nruns);
  }
  return alpha;
}


/** Gaussian neighbourhood (std shrinks exponentially with t).
	 
 */
double som_neighbourhood_gaussian( int x, int bmu, struct som_struct *s, int t ){
  double sigma;	  
  double h;

  /* sigma shrinks with runs, fast during initial phase, then slower */
  sigma = exp(-(double)t/(double)s->initial_runs)*(s->n/2) + 1;
  
  h =  exp(-pow( s->connectivity[x][bmu], 2.0)/(double)(2.0*pow(sigma, 2.0)));

  return h;
}

/** apply SOM-training rule
	 \f[
	 m_i(t+1) = m_i(t) + h_{ci}(t)[ x(t) - m_i(t) ]
	 \f]
	 
 */
void som_train_from_data( Som *s, double **X, int dim, int nsamples ){
  int i, j, t;
  int bmu=0; /* best matching unit */
  double h, tmp; 
  double bmu_score=0; /* best matching unit score */
  double *input;
  double alpha;

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
		alpha = s->time_decay( t,s->nruns,s->initial_runs );
		for( j=0; j<dim; j++ ){
		  s->m[i][j] += alpha*h*(input[j]-s->m[i][j]);
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
