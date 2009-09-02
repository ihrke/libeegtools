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
 * \brief Read EEG-data from files.

 * \todo implement reader for EEGlab set-files
 */
#ifndef READER_H
# define READER_H
#include "definitions.h"
#include "mathadd.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <strings.h>

/*typedef unsigned char uint8;*/

#define MAX_LINE_LENGTH 500

/** Matlab Defintions to read .mat files */
#define miINT8    1 /* 8 bit, signed */
#define miUINT8   2 /* 8 bit, unsigned */
#define miINT16   3 /* 16-bit, signed */
#define miUINT16  4 /* 16-bit, unsigned */
#define miINT32   5 /* 32-bit, signed */
#define miUINT32  6 /* 32-bit, unsigned */
#define miSINGLE  7 /* IEEE 754 single format */
#define miDOUBLE  9 /* IEEE 754 double format */
#define miINT64  12 /* 64-bit, signed */
#define miUINT64 13 /* 64-bit, unsigned */
#define miMATRIX 14 /* MATLAB array */

#define mxCELL_CLASS   1  /* Matlab Cell-Array */
#define mxSTRUCT_CLASS 2  /* Matlab Struct */
#define mxOBJECT_CLASS 3  /* Matlab Object */
#define mxCHAR_CLASS   4  /* Matlab Object */
#define mxSPARSE_CLASS 5  /* Matlab Object */
#define mxDOUBLE_CLASS 6  /* Matlab Object */
#define mxSINGLE_CLASS 7  /* Matlab Object */
#define mxINT8_CLASS   8  /* Matlab Object */
#define mxUINT8_CLASS  9  /* Matlab Object */
#define mxINT16_CLASS  10 /* Matlab Object */
#define mxUINT16_CLASS 11 /* Matlab Object */
#define mxINT32_CLASS  12 /* Matlab Object */
#define mxUINT32_CLASS 13 /* Matlab Object */

#ifdef __cplusplus
extern "C" {
#endif
  /** \weakgroup reader
		\{ */
  /* -------------------------------------------------------------------- */
  /* functions */
  /* -------------------------------------------------------------------- */
  EEGdata*        read_continuous_eeg_from_binfile(const char *file, int C, int n);
  double**        read_double_matrix_ascii(const char *fname, int xdim, int ydim, double **d);
  EEGdata_trials* read_eegtrials_from_raw(const char *file);
  double*         read_double_vector_ascii( const char *fname, int N, double *v );


#ifdef MATIO
  /* eeglab/matlab file reader functions if matio library is installed */
  EEG* read_eeglab( FILE *f );  
  EEG* read_eeglab_file( const char *file );
#endif


  char* read_line( FILE *f, char *line );
  ChannelInfo* read_chaninfo_ced( const char *fname, ChannelInfo *chans );

  /* matlab readers not functional, yet */
  EEGdata_trials* read_segmented_eeg_from_eeglabset(const char *file);
  int is_compressed_format(uint32_t first, int swapflag);
  /** \} */
#ifdef __cplusplus
}
#endif

#endif
