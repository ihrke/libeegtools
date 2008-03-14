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

/**\file averaging.h
 * Averaging functions.
 * \defgroup averaging Averaging functions
 *\{
 *    \defgroup timewarp Timewarping-functions
 *\}
 */
#ifndef AVERAGING_H
# define AVERAGING_H
#include "mathadd.h"
#include "definitions.h"

/* ---------------------------------------------------------------------------- 
   -- Timewarping                                                            -- 
   ---------------------------------------------------------------------------- */
/** \addtogroup timewarp
 *\ingroup averaging
 *\{
 */
double* ADTW (const double *s1, int n1, const double *s2, int n2, double *avg);
double* PADTW(const double **s, int N, int n, int zero, int *sR, double *wa);

double timewarp(const double *ref, int n1, double *d, int n2, 
		double theta1, double theta2);
double* avgwarp_from_path(const double *u, int J, const double *s, int K, 
			  const WarpPath *P, double *avg);

double get_warppath(const double *u, int J, const double *s, int K,
		    double theta1, double theta2, int* path);
WarpPath* get_warppath2(const double *u, int J, const double *s, int K,
			double theta1, double theta2);
WarpPath* get_warppath3(const double *u, int J, 
			const double *s, int K,	
			double theta1, double theta2, double *Djk);

void   warp_to_path(double *s, int K, int *path, int J);
double* warpavg_two(const double *u, const double *s, int n, int zero, 
		    int sRu, int sRs, double *wavg);
double* warpaverage(const double **ui, int N, int n,
		    const int **markers, double *wa);

double* iterative_warpavg(const ModelData *d, double *hatu);
double** diffmatrix(ModelData *m, double **dm);
/** \} */

#endif
