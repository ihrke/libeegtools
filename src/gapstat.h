/***************************************************************************
 *   Copyright (C) 2008 by Matthias Ihrke                                  *
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


/**\file gapstat.h
 * \brief \ref status_inprogress Gap-statistic for determining a suitable
        number of clusters from the data.

	 \experimental_do_not_document 

	  Reference paper: 
	  \verbatim
	  @article{tibshirani2002estimating,
  	    title={{Estimating the number of clusters in a data set via the gap statistic}},
		 author={Tibshirani, R. and Walther, G. and Hastie, T.},
		 journal={Journal of the Royal Statistical Society: Series B (Statistical Methodology)},
		 volume={63},
		 number={2},
		 pages={411--423},
		 year={2002},
		 publisher={John Wiley \& Sons}
	 }
	 \endverbatim

	  \todo rewrite with Array struct

 */

#ifndef GAPSTAT_H
# define GAPSTAT_H
#include <time.h>
#include "mathadd.h"
#include "definitions.h"
#include "array.h"
#include "averaging.h"
#include "warping.h"
#include "distances.h"
#include "helper.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct {
	 int K;            /**< maximum number of clusters */
	 int B;            /**< Monte-Carlo repetition of reference distr. calc. */
	 double *gapdistr; /**< final Gap-Statistic Distribution */
	 double *sk;       /**< modified standard deviation for gapdistr */
	 double *Wk;       /**< within scatter for data */
	 double **Wkref;   /**< within scatter for reference-dist */
	 int khat;         /**> best number of clusters */
	 ProgressBarFunction progress; /**> progress-bar callback */
  } GapStatistic;


#if 0

  GapStatistic* gapstat_init( GapStatistic *g, int K, int B );
  void          gapstat_free( GapStatistic *g );
  void          gapstat_print( FILE *out, GapStatistic *g );
  void          gapstat_calculate( GapStatistic *gap, double **X, int n, int p, 
											  VectorDistanceFunction distfunction, const double** D );
  
  double** gap_get_reference_distribution_simple( const double **X, int n, int p, double **Xr );
  double** gap_get_reference_distribution_svd   ( const double **X, int n, int p, double **Xr );

  int      eeg_best_num_clusters_gapstat( const EEG *eeg, VectorDistanceFunction distfunction,
														OptArgList *optargs );
#endif

#ifdef __cplusplus
}
#endif


#endif
