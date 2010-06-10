
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
#include "chaninfo.h"
#include "array.h"
#include "slist.h"

#ifdef __cplusplus
extern "C" {
#endif

  /** \brief describe a type of time-marker.
	*/
  typedef struct{
	 uint class; /**< used to group different types of markers (e.g. different stimuli) */
	 char label[MAX_LABEL_LENGTH]; /**< a label for the marker */
  } MarkerType;

  /** \brief a complete EEG-data set with meta-inforation.
	*/
  typedef struct{
	 char         *filename;
	 char         *comment;
	 uint         nbchan;  /**< number of channels */
	 uint         ntrials; /**< number of trials = dim(eeg) */ 
	 uint         n;       /**< number of samples */
    double       sampling_rate; /**< in Hz */

	 Array       *times;  /**< times vector (1D DOUBLE Array); n-long */
	 ChannelInfo  *chaninfo; /**< location and other information
										 about the channels */

	 Array         *data;   /**< 3D DOUBLE Array: channels x trials x samples */

	 Array         *markers; /**< 2D UINT Array: trials x max_num_markers; 
									      this field gives the offset in sampling points
									      of each time-marker (e.g. stimulus onset,...).
									 */
	 Array         *marker_type; /**< 2D UINT Array: trials x max_num_markers;
											  this field describes the type of the marker, by
										     indexing into the EEG->marker_types field that
										     contains the complete description of the 
										     available markers
										  */
	 SingleList    *marker_types; /**< contains the complete description of the 
												available markers; the list->content field 
											   is a \ref MarkerType struct */
	 void *additional;            /**< arbitrary, additional information for the
												data set. This can be anything e.g. ICA-components,
												distance matrix or even both. The additional data
												is automatically saved/retrieved from files (if 
												supported by the format). To use it in your program
												you wil of course need to know what it is.
											*/
	 ulong nbytes_additional;     /**< number of bytes in EEG->additional */
  } EEG;

  
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
  char* eeg_sprint( char *out, const EEG *eeg, int preview );
  void eeg_print( FILE *out,  const EEG *eeg, int preview );
  EEG* eeg_clone( const EEG *eeg, int flags );
  void eeg_free ( EEG *eeg );

  void print_channelinfo( FILE* out, const ChannelInfo *c );

#ifdef __cplusplus
}
#endif

#endif
