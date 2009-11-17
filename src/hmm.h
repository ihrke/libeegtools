/***************************************************************************
 *   Copyright (C) 2008/2009 by Matthias Ihrke                                  *
 *   mihrke@uni-goettingen.de                                              *
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

/**\file hmm.h
 \brief Hidden Markov Models.

 \warning For now this is only experimental but plans are to complete this file.

 Implements the continuous profile model from:

 Listgarten et al. Multiple alignment of continuous time
 series. Advances in Neural Information Processing Systems (2005)
 vol. 17 pp. 817–824

 */
#ifndef HMM_H
# define HMM_H

#include "mathadd.h"
#include "definitions.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**\ingroup cphmm
	*\{
	*/ 

  /** Continuous Profile Model struct. See Listgarten et al. 2005.
	*/
  typedef struct{
	 int K; /**< number of trials */
	 int n; /**< number of samples per trial */
	 int M; /**< number of samples in latent trace */

	 int Q; /**< number of amplitude states */
	 double *q; /**< (q1,...,q_Q) amplitude factor states */
  
	 double *u; /**< K-vector; global amplitude multiplication factor */
  
	 int J; /**< maximum number of jumps in warping*/
	 double *d; /**< (d1, ..., di, ..., dJ); probabilities to do a jump of length i */
	 double s[2]; /**< (s0, s1); probabilities to go up/down one amplitude state */

	 double **tau; /**< (K x N) warping for each trial w.r.t. z */
	 int    **phi; /**< (K x N) amplitude states for each trial; q(phi^k_i) is the actual amplitude */

	 double *z; /**< M-vector; latent trace; */
	 double sigma; /**< noise factor; standard deviation in emission prob. density */

	 gsl_rng_type *random_number_type; /**< GSL-random number generator type */
	 gsl_rng *rng; /**< the GSL random number generator */
	 long seed;   /**< seed for the rng */
  } CPHiddenMarkovModel;

  double cphmm_get_transition_prob( CPHiddenMarkovModel *m, int k, int i, int j );


  CPHiddenMarkovModel* cphmm_alloc( int K, int n, int M, int Q, int J );
  CPHiddenMarkovModel* eegtrials_cphmm_init( EEG *eeg, int channel, double **X );
  void cphmm_init( CPHiddenMarkovModel *m, double **X );
  void cphmm_free( CPHiddenMarkovModel *m );

  /*\}*/

#ifdef __cplusplus
}
#endif


#endif
