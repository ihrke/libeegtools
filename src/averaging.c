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

/**\file 
	\copydoc averaging.h
 */
#include "averaging.h"
#include "clustering.h"
#include "eeg.h"
#include "linalg.h"

/** \brief example function to pass as SignalAverageFunction doing a pointwise average.

	 This function is only for testing purposes. The "hierarchically" computed
	 pointwise average should be the same as the normal average.

	 \param input data (N x C x n)
	 \param index average data[idx[0]] and data[idx[1]]
	 \param weights number of trials "behind" each average
	 \param optional arguments
	 \return the average
*/
Array* average_example( const Array *data, uint idx[2], double weights[2], OptArgList *optargs ){
  int i,j;
  if( data->ndim!=3 || data->dtype!=DOUBLE ){
	 errprintf("data must be 3D and DOUBLE, have %i, %i\n", data->ndim, data->dtype );
	 return NULL;
  } 
  int N,C,n;
  N = data->size[0];
  C = data->size[1];
  n = data->size[2];
  Array *avg = array_new2( DOUBLE, 2, C, n );
  for( i=0; i<C; i++ ){
	 for( j=0; j<n; j++ ){
		mat_IDX( avg, i, j )=
		  (weights[0]*array_INDEX3( data, double, idx[0], i, j ))+
		  (weights[1]*array_INDEX3( data, double, idx[1], i, j ));
	 }
  }									

  return avg;
}
Array* average_warpmarkers( const Array *data, uint idx[2], double weights[2], OptArgList *optargs ){
}

/** \brief Calculate a hierarchical average based on a cluster-analysis.
	 
	 The data is average based on the first dimension of the data.
	 I.e., the output is C x n dimensional.

	 A cluster analysis is carried out using the distance-matrix and the 
	 resulting dendrogram is followed to average two trials at a time:

	 \image html dendrogram.jpg

	 The averaging uses a callback-function to calculate the average of two signals
	 at a time. This function must be passed by the use of a SignalAverageFunction .

	 \param data 2D or 3D double-array containing: 
	      - N trials (dimension 1)
			- of dimensionality C (dimension 2); e.g. channels
			- and n time-points (dimension 3)
			if the data is 2D, we assume C=1			
	\param distmat a distance matrix (N x N) given distances between trials in data
	\param avgfct a callback-function to calculate the average of two Cxn signals at a time
	\param optargs may contain:
  	 - <b>Directly used by this function:</b>
	 - <tt>linkage=void*</tt> rule for building the dendrogram, default=\c dgram_dist_completelinkage
	 - <tt>progress=void*</tt> progress-bar function, called every now and then.
	 - <b>You might want to pass additional optargs for the SignalAverageFunction and the 
	   linkage function</b>
	 \return pointer to final average (Cxn)
 */
Array* hierarchical_average( const Array *data, const Array *distmat, 
									  SignalAverageFunction avgfct, OptArgList *optargs ){
#if 1
  Dendrogram *T, *Tsub; 
  ProgressBarFunction progress=NULL;  
  LinkageFunction linkage=dgram_dist_completelinkage;
  void *ptr;

  bool ismatrix;
  matrix_CHECKSQR( ismatrix, distmat );
  if( !ismatrix ) return NULL;
  if( data->ndim>3 || data->ndim<2 || data->dtype!=DOUBLE ){
	 errprintf( "Data needs to be 2D or 3D, got %iD and of type DOUBLE\n", data->ndim );
	 return NULL;
  }

  optarg_PARSE_PTR( optargs, "linkage", linkage, LinkageFunction, ptr );
  optarg_PARSE_PTR( optargs, "progress", progress, ProgressBarFunction, ptr );

  int N,C,n;
  N = data->size[0];
  C = (data->ndim==2)?1:(data->size[1]);
  n = (data->ndim==2)?(data->size[1]):(data->size[2]);
  dprintf("Data is (%i x %i x %i)\n", N, C, n );

  /* wrapper for 2D-data -> work on d instead of data */
  Array *d=array_copy( data, TRUE );
  if( d->ndim==2 ){
	 d->ndim=3;
	 free( d->size );
	 d->size=(uint*)malloc( 3*sizeof(uint) );
	 d->size[0]=N; d->size[1]=C; d->size[2]=n;
  }

  /*------------------- computation --------------------------*/
  /* build hierarchical dendrogram */
  T = agglomerative_clustering( distmat, linkage );

  int i;
  double weights[2];
  int *indices = (int*)malloc( N*sizeof(int) );
  for( i=0; i<N; i++ ){
	 indices[i] = 1;
  }  

  /* now walk the tree to find pairs of trials to match */
  Array *new, 
	 *avg=array_new2( DOUBLE, 2, C, n );
  int trials_left = N;
  uint idx[2]={0,0};
  if( progress ) progress( PROGRESSBAR_INIT, N );

  while( trials_left >= 2 ){
	 if( progress ) progress( PROGRESSBAR_CONTINUE_LONG, N-trials_left );
	 
	 Tsub = dgram_get_deepest( T );
	 idx[0] = Tsub->left->val;
	 idx[1] = Tsub->right->val;
	 
	 weights[0] = indices[idx[0]]/(double)(indices[idx[0]]+indices[idx[1]]);
	 weights[1] = indices[idx[1]]/(double)(indices[idx[0]]+indices[idx[1]]);
	 /* Array* average_example( const AADTWrray *data, uint idx[2], 
		       double weights[2], OptArgList *optargs ) */
	 new=avgfct( d, idx, weights, optargs );
	 memcpy( array_INDEXMEM3( d, idx[0], 0, 0 ),
				new->data, C*n*sizeof( double ) );
	 array_free( new );

	 indices[idx[0]] += indices[idx[1]];
	 indices[idx[1]]=-1; 			  /* never to be used again */
	 
	 /* replace node by leaf representing average(idx1, idx2) */
	 Tsub->val = idx[0];
	 Tsub->left = NULL;
	 Tsub->right = NULL;	
	 trials_left--;
  }
  /* idx[0] is the final average */
  memcpy( avg->data, array_INDEXMEM3( d, idx[0], 0, 0 ),
			 C*n*sizeof( double ) );
  /*------------------- /computation --------------------------*/

  free( indices );
  array_free( d );

  array_dimred( avg );
  return avg;
#endif
}



/** Calculate a simple, pointwise average across trials
	 \f[
	 \hat{s}(t) = \frac{1}{N}\sum_{i=1}^{N}s_i(t)
	 \f]
	 \param eeg input
	 \return freshly allocate EEG-struct
 */
EEG*     eeg_simple_average( const EEG *eeg ){
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
	 /* TODO */
  }
  return out;
}
/** Calculate a pointwise average across channels
	 \f[
	 \hat{s}_i(t) = \frac{1}{C}\sum_{c=1}^{C}s^c_i(t)
	 \f]
	 \param eeg input
	 \return freshly allocate EEG-struct
 */
EEG*     eeg_average_channels( const EEG *eeg ){
  EEG *out;
  int c, i, j;

  out = eeg_init( 1, eeg->ntrials, eeg->n ); /* the average */
  out->sampling_rate=eeg->sampling_rate;
  if( eeg->times ){
	 out->times = (double*) malloc( eeg->n*sizeof(double) );
	 memcpy( out->times, eeg->times, eeg->n*sizeof(double) );
  }
  eeg_append_comment( out, "output from eegtrials_average_channels()\n" );

  for( i=0; i<eeg->ntrials; i++ ){
	 for( j=0; j<eeg->n; j++ ){
		out->data[0][i][j] = 0;
		for( c=0; c<eeg->nbchan; c++ ){
		  out->data[0][i][j] += eeg->data[c][i][j];
		}
		out->data[0][i][j] /= (double)eeg->nbchan;
	 }
  }
  return out;
}
