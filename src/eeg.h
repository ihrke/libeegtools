
/* **************************************************************************
 *   Copyright (C) 2009 by Matthias Ihrke   *
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

/**\file eeg.h
 * \brief \ref status_stable EEG struct handling functions
 */
#ifndef EEG_H
# define EEG_H

#include "definitions.h"
#include <string.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

#define EEG_CLONE_ALL        0     /**< clone everything */
#define EEG_CLONE_NODATA     2<<1  /**< clone everything except 
												  the eeg->data field */
#define EEG_CLONE_NOMARKERS  2<<2  /**< clone everything except
												  the eeg->markers, eeg->nmarkers
												  eeg->marker_labels fields */
#define EEG_CLONE_NOCHANINFO 2<<3  /**< clone everything except
												  the eeg->chaninfo field */

  
  EEG* eeg_init            ( int nbchan, int ntrials, int nsamples );
  EEG* eeg_init_markers    ( int nmarkers_per_trial, EEG *eeg );
  EEG* eeg_extract_channels( EEG* eeg, const int *channels, int nchannels, bool alloc );
  EEG* eeg_extract_trials  ( EEG* eeg, const int *trials,   int ntrials,   bool alloc );

  void eeg_append_comment( EEG *eeg, const char *comment );
  void eeg_print( FILE *out,  const EEG *eeg, int preview );
  EEG* eeg_clone( const EEG *eeg, int flags );
  void eeg_free ( EEG *eeg );

  void print_channelinfo( FILE* out, const ChannelInfo *c );

#ifdef __cplusplus
}
#endif

#endif
