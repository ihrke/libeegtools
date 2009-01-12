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
 *    \defgroup otheravg Other Averaging functions
 *\}
 */
#ifndef AVERAGING_H
#define AVERAGING_H
#include "mathadd.h"
#include "definitions.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \addtogroup otheravg
 *\ingroup averaging
 *\{
 */
  double* simple_average_2v(const double *s1, const double *s2, int n, double *avg);
  double* simple_average_nv(const double **s, int N, int n, double *avg);
  double* alternate_average_nv(const double **s, int N, int n, double *avg);

  EEGdata* eegtrials_simple_average( EEGdata_trials *eeg, EEGdata *avg);
  double*  eegtrials_simple_average_channel( EEGdata_trials *eeg, int channel, double *avg );

/** \} */
#ifdef __cplusplus
}
#endif

#endif
