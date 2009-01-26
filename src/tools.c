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

/** Remove baseline from eeg-data. Computes the average of the eeg-values 
	 in the given time-window and substracts it from the all data.
	 \param eeg - the struct is overwritten
	 \param times - eeg->n long; indicates which values in eeg->d corresponds 
	                a which time
	 \param win_from - the time-window in ms (starting from)
	 \param win_to - the time-window in ms (until)
	 \note both win_from and win_to must be in the times-array
*/
void eeg_remove_baseline( EEGdata *eeg, const double *times, 
								  double win_from, double win_to ){
  int lim[2], chan;
  double mean;

  lim[0] = closest_index( times, eeg->n, win_from );
  lim[1] = closest_index( times, eeg->n, win_to );
  if( lim[1]<=lim[0] ){
	 errprintf( "The specified limits are funny: %.2f-%.2f (results in %i-%i)\n", 
					win_from, win_to, lim[0], lim[1] );
  }
  for( chan=0; chan<eeg->nbchan; chan++ ){
	 mean = gsl_stats_mean( eeg->d[chan]+lim[0], 1, lim[1]-lim[0] );
	 // printf(" mean=%f\n", mean);
	 vector_minus_scalar( eeg->d[chan], eeg->n, mean );
  }
}

/** Remove baseline from eeg-data. Computes the average of the eeg-values 
	 in the given time-window and substracts it from the all data. It's done
	 for all trials in the struct calling  eeg_remove_baseline().
	 \param win_from - the time-window in ms (starting from)
	 \param win_to - the time-window in ms (until)
	 \note both win_from and win_to must be in the times-array of eeg
*/
void eegtrials_remove_baseline( EEGdata_trials *eeg, double win_from, double win_to ){
  int t;
  for( t=0; t<eeg->ntrials; t++ ){
	 dprintf(" Baseline for trial %i\n", t );
	 eeg_remove_baseline( eeg->data[t], eeg->times, win_from, win_to );
  }
}
