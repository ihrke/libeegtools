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

/**\file clustering.h
 * Clustering functions.
 * Groups:
 \defgroup clustering Clustering functions
 \{
    \defgroup distances Distance measures
    \defgroup gap Gap-statistic
 \}
 */
#ifndef CUSTERING_H
# define CLUSTERING_H
#include <time.h>
#include "mathadd.h"
#include "definitions.h"
#include "averaging.h"

#ifdef __cplusplus
extern "C" {
#endif

/**\addtogroup clustering
 *\{
 */
double** diffmatrix(ModelData *m, double **dm);
double** eegtrials_diffmatrix_channel(EEGdata_trials *eeg, 
												  double(*dist)(EEGdata*,EEGdata*,int), 
												  int channel);

Clusters* kmedoids(const double **dist, int N, int K);
void free_cluster(Clusters *c);
void print_cluster(Clusters *c);
Clusters* init_cluster(int K, int maxN);
void copy_cluster(Clusters *dest, const Clusters *src);
void diffmatrix_standardize(double **d, int N);
/** \} */

/**\addtogroup distances 
	\ingroup clustering
	\{*/
double clusterdist_euclidean_pointwise(EEGdata  *s1, EEGdata *s2, int channel);
double clusterdist_tw_complete(EEGdata *s1, EEGdata *s2, int channel);
double clusterdist_tw_markers(EEGdata *s1, EEGdata *s2, int channel);
/** \} */


/**\addtogroup gap
	\ingroup clustering
	\{
	Reference paper: Tibshirani, 2001
*/
int      gap_get_K(const double *gapstat, int k);
double*  gap_get_gapstat(const double *Wk, const double **Wkref, int B, int k, double *gapstat);
double** gap_get_reference_distribution(const double **d, int N, int n, double **ref);
double   gap_get_within_scatter(const double **d, int N, const Clusters *c);
double*  gap_get_within_scatter_distribution(const double **d, int N, int k, const Clusters **c, double *Wk);
/** \} */


#ifdef __cplusplus
}
#endif


#endif
