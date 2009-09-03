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
 * \brief EEG struct handling functions
 */
#ifndef EEG_H
# define EEG_H

#include "definitions.h"
#include <string.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif


  EEG* eeg_init            ( int nbchan, int ntrials, int nsamples );
  EEG* eeg_extract_channels( EEG* eeg, int *channels, int nchannels, Boolean alloc );
  EEG* eeg_extract_trials  ( EEG* eeg, int *trials,   int ntrials,   Boolean alloc );

  void eeg_print( FILE *out,  EEG *eeg, int preview );
  EEG* eeg_clone( EEG *eeg );

 /**\addtogroup helperstruct
  *\ingroup helper
  *\{*/
  /* ---------------------------------------------------------------------------- 
	  -- Helper functions (struct handling)                                     -- 
	  ---------------------------------------------------------------------------- */
  /* constructors */
  EEGdata*        init_eegdata(int nbchan, int nsamples, int nmarkers);
  void            reset_eegdata( EEGdata* eeg );
  EEGdata_trials* init_eegdata_trials(int nbtrials, int markers_per_trial, 
												  int nbchan, int nbsamples, double *times);
 
  
  /* destructors */
  void    free_eegdata(EEGdata *eeg);
  void    free_eegdata_trials(EEGdata_trials *eeg);

  /* convenience functions for structs */
  int       eegdata_cmp_settings( const EEGdata *s1, const EEGdata *s2 );
  int       copy_similar_eegdata( EEGdata *dest, const EEGdata *source );
  EEGdata*         clone_eegdata( const EEGdata *source );
  EEGdata_trials*  clone_eegdata_trials( const EEGdata_trials *source );

  /* printing-functions */
  void    print_eegdata_trials(FILE *out, const EEGdata_trials *eeg);
  void    print_eegdata(FILE *out, const EEGdata *eeg);
  void    print_channelinfo( FILE *out, const ChannelInfo *c );
  /**\}*/


#ifdef __cplusplus
}
#endif

#endif
