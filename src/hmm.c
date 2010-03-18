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

#include "hmm.h"


/** \ref status_inprogress
	 return state transition probabilities from state i to state j in trace k:
	 \f[
	 p^k(\pi_i|\pi_{i-1}) = p^k(\phi_i|\phi_{i-1})p^k(\tau_i|\tau_{i-1})
	 \f]
	 where
	 \f[
	 \pi_i = (\phi_i, \tau_i)
	 \f]
	 Details see Listgarten et al. 2005
	 \todo not implemented yes
 */
double cphmm_get_transition_prob( CPHiddenMarkovModel *m, int k, int i, int j ){
  double p_phi=0, p_tau=0;


  if( j-i >= 0 && j-i < m->J ){
	 p_tau = m->d[j-i];
  }
  
  return p_phi*p_tau;
}


/** allocate memory for the CPM.
	 \param K number of trials
	 \param n  number of samples per trial 
	 \param M number of samples in latent trace 
	 \param Q number of amplitude states
	 \param J maximum number of jumps in warping
 */
CPHiddenMarkovModel* cphmm_alloc( int K, int n, int M, int Q, int J ){
  CPHiddenMarkovModel *m;
  int i;

  m = (CPHiddenMarkovModel*)malloc( sizeof(CPHiddenMarkovModel) );
  m->K = K;
  m->n = n;
  m->M = M;
  m->Q = Q;
  m->J = J; 
  
  m->q = (double*) malloc( Q*sizeof(double) );
  m->u = (double*) malloc( K*sizeof(double) );
  m->d = (double*) malloc( J*sizeof(double) );
  m->z = (double*) malloc( M*sizeof(double) );

  m->tau = (double**) malloc( K*sizeof(double*) );
  m->phi = (int   **) malloc( K*sizeof(int   *) );
  for( i=0; i<K; i++ ){
	 m->tau[i] = (double*) malloc( n*sizeof(double) );
	 m->phi[i] = (int   *) malloc( n*sizeof(int   ) );
  }

  gsl_rng_env_setup();
  m->random_number_type = (gsl_rng_type *)gsl_rng_default;
  m->rng = gsl_rng_alloc (m->random_number_type);

  return m;
}

/** initializing the model according to the parameters from
	 Listgarten et al. 2005:
	 
	 - sigma = 15% of max(X[0])-min(X[0])
	 - z = first observed time-series + gaussian noise (0, sigma)
	       upsampling to match M
	 - u = (1,...,1)
	 - d,s = uniform
	 - q is filled with Q numbers evenly spaced in log-space [0,2]
	 \param m the model
	 \param X the data, K x n
 */
void cphmm_init( CPHiddenMarkovModel *m, double **X ){
  int i;
  double minx, maxx;
  int start;

  minx = dblp_min( X[0], m->n, NULL );
  maxx = dblp_max( X[0], m->n, NULL );
  m->sigma = 0.15*(maxx-minx);

  start = (m->M-(2*m->n))/2.0;
  if( start<0 ){
	 warnprintf("M<2*n: %i<%i -> skipping z-initialization\n", m->M, 2*m->n);
  }

  for( i=0; i<start; i++ ){
	 m->z[i] = minx;
  }
  for( i=start; i<2*m->n; i++ ){
	 m->z[i] = X[0][i-start]+gsl_ran_gaussian( m->rng, m->sigma );
  }

  for( i=0; i<m->K; i++ ){
	 m->u[i] = 1;
  }
  for( i=0; i<m->J; i++ ){
	 m->d[i] = i/(double)m->J;
  }
  m->s[0] = 0.5;
  m->s[1] = 0.25;

  for( i=0; i<m->Q; i++ ){
	 m->q[i] = i*(2.0/(double)m->Q);
  }

}

/** init the model from EEG-data. The data matrix X is a return value.
	 \param eeg
	 \param channel - which channel is used?
	 \param X - the converted data. X should contain K empty pointers. These are set
	            to the appropriate data in eeg (not copied).

 */
CPHiddenMarkovModel* eeg_cphmm_init( EEG *eeg, int channel, double **X ){
  CPHiddenMarkovModel *m;
  int M;
  double epsilon = 0.2;
  int Q;
  int J;
  int i;

  M = (2+epsilon)*eeg->n;
  Q = 7;
  J = 3;

  m = cphmm_alloc( eeg->ntrials, eeg->n, M, Q, J );
  for( i=0; i<m->K; i++ ){
	 X[i] = eeg->data[channel][i];
  }
  cphmm_init( m, X );

  return m;
}

/** deallocate the model's struct. RNG is dealloc'ed as well.
 */
void cphmm_free( CPHiddenMarkovModel *m ){
  int i;
  if(m->q) free(m->q);
  if(m->u) free(m->u);
  if(m->d) free(m->d);
  if(m->z) free(m->z);
  for( i=0; i<m->K; i++ ){
	 if( m->tau[i] )
		free( m->tau[i] );
	 if( m->phi[i] )
		free( m->phi[i] );
  }
  if( m->tau ) 
	 free( m->tau );
  if( m->phi )
	 free( m->phi );

  gsl_rng_free (m->rng);
}
