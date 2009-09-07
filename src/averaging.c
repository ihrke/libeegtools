/***************************************************************************
 *   Copyright (C) 2008 by Matthias Ihrke   *
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

#include "averaging.h"
#include "eeg.h"

/** Calculate a simple, pointwise average across trials
	 \f[
	 \hat{s}(t) = \frac{1}{N}\sum_{i=1}^{N}s_i(t)
	 \f]
	 \param eeg input
	 \return freshly allocate EEG-struct
 */
EEG*     eeg_simple_average( EEG *eeg ){
  EEG *out;
  int c;

  out = eeg_init( eeg->nbchan, 1, eeg->n ); /* the average */
  out->sampling_rate=eeg->sampling_rate;
  if( eeg->times ){
	 out->times = (double*) malloc( eeg->n*sizeof(double) );
	 memcpy( out->times, eeg->times, eeg->n*sizeof(double) );
  }
  if( eeg->chaninfo ){
	 out->chaninfo = (ChannelInfo*) malloc( eeg->nbchan*sizeof(ChannelInfo) );
	 memcpy( out->chaninfo, eeg->chaninfo,  eeg->nbchan*sizeof(ChannelInfo) );
  }
  eeg_append_comment( out, "output from eegtrials_simple_average()\n" );
  for( c=0; c<eeg->nbchan; c++ ){
	 out->data[c][0] = simple_average_nv( eeg->data[c], eeg->ntrials, eeg->n, out->data[c][0] );
  }
  return out;
}

/** Calculate alternate average of n vectors
	 \f[
	 \hat{s}(t) = \frac{1}{N}\sum_{i=1}^{N} (-1)^{i-1}s_i(t)
	 \f]
	 \param eeg input
	 \return freshly allocate EEG-struct
 */
EEG*     eeg_alternate_average( EEG *eeg ){
  EEG *out;
  int c;

  out = eeg_init( eeg->nbchan, 1, eeg->n ); /* the average */
  eeg_append_comment( out, "output from eegtrials_simple_average()\n" );
  for( c=0; c<eeg->nbchan; c++ ){
	 out->data[c][0] = alternate_average_nv( eeg->data[c], eeg->ntrials, eeg->n, out->data[c][0] );
  }
  return out;
}


/** Calculate a simple, pointwise average of two vectors. 
	 \f[ \hat{s}(t) = \frac{(s_1(t)+s_2(t))}{2} \f]
	 \param avg - user-allocated memory of length n; if NULL, the function 
	              allocates the memory
 */
double* simple_average_2v(const double *s1, const double *s2, int n, double *avg){
  int i;
  if(!avg)
    avg = (double*)malloc(n*sizeof(double));
  for(i=0; i<n; i++)
    avg[i] = (s1[i]+s2[i])/2.0;
  return avg;
}

/** Calculate a simple, pointwise average of n vectors
	 \f[
	 \hat{s}(t) = \frac{1}{N}\sum_{i=1}^{N}s_i(t)
	 \f]
	 \param s - data N x n (trials x samples)
	 \param N trials
	 \param n samples
 * \param avg - user-allocated memory of length n; if NULL, the function 
 *              allocates the memory
 */
double* simple_average_nv(const double **s, int N, int n, double *avg){
  int i, j;
  if(!avg)
    avg = (double*)malloc(n*sizeof(double));
  for(i=0; i<n; i++) {
	 avg[i] = 0;
	 for( j=0; j<N; j++ ){
		avg[i] += s[j][i];
	 }
	 avg[i] /= (double) N;
  }
  return avg;
}

/** Calculate alternate average of n vectors
	 \f[
	 \hat{s}(t) = \frac{1}{N}\sum_{i=1}^{N} (-1)^{i-1}s_i(t)
	 \f]
 * \param avg - user-allocated memory of length n; if NULL, the function 
 *              allocates the memory
 */
double* alternate_average_nv(const double **s, int N, int n, double *avg){
  int i, j;
  if(!avg)
    avg = (double*)malloc(n*sizeof(double));
  for(i=0; i<n; i++) {
	 avg[i] = 0;
	 for( j=0; j<N; j++ ){
		avg[i] += pow(-1, j-1)*s[i][j];
	 }
	 avg[i] /= (double) N;
  }
  return avg;
}
