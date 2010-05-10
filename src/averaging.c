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

Array* average_example( const Array *data, uint idx[2], double weights[2], OptArgList *optargs ){
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
#if 0
  Dendrogram *T, *Tsub; 
  ProgressBarFunction progress=NULL;  
  LinkageFunction linkage=dgram_dist_completelinkage;
  void *ptr;
  if( optarglist_has_key( optargs, "linkage" ) ){
	 ptr = optarglist_ptr_by_key( optargs, "linkage" );
	 if( ptr ) linkage = (LinkageFunction)ptr;
  }
  /* build hierarchical dendrogram */
  T = agglomerative_clustering( (const double**)distmatrix, eeg->ntrials, linkage );
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
