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

#ifndef MATH_ADD_H
# define MATH_ADD_H
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include "helper.h"
#include "definitions.h"

/** \ingroup helpermath*/
#define MAX(a,b) ((a) > (b) ? (a):(b))
/** \ingroup helpermath*/
#define MIN(a,b) ((a) < (b) ? (a):(b))
/** \ingroup helpermath*/
#define SQR(a) ((a)*(a))


/* -------- Math ----------- */
/**\addtogroup helpermath
 *\{*/
double  glog(double v, int b);  
double  mad(const double *data, int n); 
int     abscmp(const void *p1, const void *p2);
double  vnorm(const double *v, int n, int p);
double  maxel(double *v, int n);
int     maxeli(int *v, int n);
int     closest_index(double *v, int n, double c);
double* sampled_line(double *ntimes, int n, double start, double end);
double* lininterp(const double *x1, const double *y1, int n1, 
						const double *x2,       double *y2, int n2);
int* linspace(int first, int last);
double* loocv(const ModelData *m, double* err,
				  double*(*model)(const ModelData*,double*));

/* ---------------------------------------------------------------------------- 
   -- Merit Measures                                                         -- 
   ---------------------------------------------------------------------------- */
double rmse(const double *r, const double *d, int n);
double snr (const double *r, const double *d, int n);

/**\}*/

#endif
