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
		avg[i] += s[i][j];
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


EEGdata* eegtrials_simple_average( EEGdata_trials *eeg, EEGdata *avg){
  int c, nbchan, n;
  
  nbchan = eeg->data[0]->nbchan;
  n = eeg->data[0]->n;
  if( avg==NULL ){
	 avg = init_eegdata( nbchan, n, 0 );
  }
  
  for( c=0; c<nbchan; c++ ){
	 eegtrials_simple_average_channel( eeg, c, avg->d[c] );
  }

  return avg;
}

double* eegtrials_simple_average_channel( EEGdata_trials *eeg, int channel, double *avg ){
  int t, i,n, N;

  n =  eeg->data[0]->n;
  N = eeg->ntrials;

  if( avg==NULL ){
	 avg = (double*)malloc( n * sizeof(double) );
  }

  memset( avg, 0, n*sizeof(double) );

  for( i=0; i<n; i++ ){
	 for( t=0; t<N; t++ ){
		avg[i] += eeg->data[t]->d[channel][i];
	 }
	 avg[i] = avg[i]/(double)N;
  }

  return avg;
}
