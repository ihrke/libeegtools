/***************************************************************************
 *   Copyright (C) 2010 by Matthias Ihrke                                  *
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

/**\file mpitools.h
 \brief \ref status_stable Tools for using LibEEGTools with MPI.
	
 */
#ifndef MPITOOLS_H
# define MPITOOLS_H
#include "mathadd.h"
#include "definitions.h"

#ifdef __cplusplus
extern "C" {
#endif

  ulonglong_t sizeof_eegstream( const EEG *eeg );
  ulonglong_t sizeof_channelinfos( const ChannelInfo *chaninfo, int nbchan );

  char* eeg_to_stream( const EEG *eeg, char *stream, ulonglong_t *n );
  EEG*  stream_to_eeg( const char *stream, ulonglong_t n, EEG *eeg );

#ifdef __cplusplus
}
#endif


#endif
