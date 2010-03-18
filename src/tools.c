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

#include "tools.h"
#include "eeg.h"
#include "mathadd.h"

/** Remove baseline from eeg-data. Computes the average of the eeg-values 
	 in the given time-window and substracts it from the all data.
	 \param eeg - the struct 
	 \param win_from - the time-window in ms (starting from)
	 \param win_to - the time-window in ms (until)
	 \note both win_from and win_to must be in the times-array of eeg->times
	 
*/
EEG* eeg_remove_baseline( EEG *eeg, double win_from, double win_to, bool alloc ){
  int lim[2], c, i;
  double mean;
  EEG *eeg_out;

  if( !eeg->times ){
	 errprintf("Need times-array, aborting...\n");
	 return NULL;
  }

  if( alloc ){
	 eeg_out = eeg_clone( eeg, EEG_CLONE_ALL );
  } else {
	 eeg_out = eeg;
  }

  lim[0] = closest_index( eeg->times, eeg->n, win_from );
  lim[1] = closest_index( eeg->times, eeg->n, win_to );
  if( lim[1]<=lim[0] ){
	 errprintf( "The specified limits are funny: %.2f-%.2f (results in %i-%i)\n", 
					win_from, win_to, lim[0], lim[1] );
  }
  for( c=0; c<eeg->nbchan; c++ ){
	 for( i=0; i<eeg->ntrials; i++ ){
		mean=dblp_mean( eeg->data[c][i]+lim[0], lim[1]-lim[0] );
		// printf(" mean=%f\n", mean);
		dblp_minus_scalar( eeg->data[c][i], eeg->n, mean );
	 }
  }

  return eeg_out;
}
