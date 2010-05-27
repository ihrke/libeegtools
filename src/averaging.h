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
 * \brief \ref status_stable Averaging functions.
 */
#ifndef AVERAGING_H
#define AVERAGING_H
#include "mathadd.h"
#include "definitions.h"
#include "array.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \weakgroup otheravg
 *\{
 */

/** \brief Function for averaging two trials from data.
	 \param input data (N x C x n)
	 \param index average data[idx[0]] and data[idx[1]]
	 \param weights number of trials "behind" each average
	 \param optional arguments
	 \return the average
 */
  typedef Array*(*SignalAverageFunction)( const Array *,uint[2],double[2],OptArgList*);
  Array* average_example( const Array *data, uint idx[2], double weights[2], 
								  OptArgList *optargs );

  Array* hierarchical_average( const Array *data, const Array *distmat, 
										 SignalAverageFunction avgfct, OptArgList *optargs );

  EEG*     eeg_simple_average   ( const EEG *eeg );
  EEG*     eeg_average_channels ( const EEG *eeg );

/** \} */
#ifdef __cplusplus
}
#endif

#endif
