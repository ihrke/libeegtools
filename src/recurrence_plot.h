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

/**\file recurrence_plot.h
 * \brief \ref status_unstable Functions using nonlinear systems-theory and recurrence-plots.
 */
#ifndef RECURRENCE_PLOT_H
# define RECURRENCE_PLOT_H
#include "mathadd.h"
#include "nonlinear.h"
#include "clustering.h"
#include "definitions.h"
#include "warping.h"
#include <gsl/gsl_sort_double.h>


/* optional args for recurrence plots */
#define RPLOT_NONE 0
#define RPLOT_FAN  1
/* #define RPLOT_DIST 2 */

#ifdef __cplusplus
extern "C" {
#endif

  /*-----------------------------------------------------------
	 - Recurrence Plots - 
	 ---------------------------------------------------------*/

  Array* recplot( const Array *s1, const Array *s2, Array *out, double epsilon,
						OptArgList *optargs );
  double* recplot_calculate_epsilons( Array *s1, Array *s2, double *eps, int fan );

  /* Line-of-synchrony detection */
  Array* recplot_los_dtw_noise( const Array *R );
  Array* recplot_los_marwan( const Array *R, int dx, int dy );

#ifdef __cplusplus
}
#endif


#endif
