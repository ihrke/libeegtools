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

/**\file reader.h
 * \brief \ref status_stable Read EEG-data from files.
 */
#ifndef READER_H
# define READER_H
#include "definitions.h"
#include "mathadd.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <strings.h>

#define MAX_LINE_LENGTH 500

#ifdef __cplusplus
extern "C" {
#endif
  /** \weakgroup reader
		\{ */
  /* -------------------------------------------------------------------- */
  /* functions */
  /* -------------------------------------------------------------------- */
  EEG*     read_continuous_eeg_from_binfile(const char *file, int C, int n);
  double** read_double_matrix_ascii(const char *fname, int xdim, int ydim, double **d);
  EEG*     read_eeg_from_raw(const char *file);
  double*  read_double_vector_ascii( const char *fname, int N, double *v );


#ifdef MATIO
  /* eeglab/matlab file reader functions if matio library is installed */
  EEG* read_eeglab_file( const char *file );
#endif

  ChannelInfo* read_chaninfo_ced( const char *fname, ChannelInfo *chans );
  /** \} */
#ifdef __cplusplus
}
#endif

#endif
