/***************************************************************************
 *   Copyright (C) 2010 by Matthias Ihrke                                  *
 *   ihrke@nld.ds.mpg.de
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
 *   aint with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/**\file io.h
 \brief \ref status_inprogress io input-output for (EEG)-data.


 These functions are for some generic data types.

 For specialized IO-functions, see
 
  - \ref io_matlab.h -- Input/Output for MATLAB (requires MatIO)
  - \ref io_wav.h -- Input/Output for WAV-Audio Files
 */
#ifndef IO_H
#define IO_H

#define MAX_LINE_LENGTH 500

#include "definitions.h"
#include "eeg.h"
#include "optarg.h"

/* include "special" IO-headers */
#include "io_matlab.h"
#include "io_wav.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* ------------- READER --------------------- */
  double** read_double_matrix_ascii(const char *fname, int xdim, int ydim, double **d);
  double*  read_dblp_ascii( const char *fname, int N, double *v );

  EEG*     read_continuous_eeg_from_binfile(const char *file, int C, int n);
  EEG*     read_eeg_from_raw(const char *file);

  ChannelInfo* read_chaninfo_ced( const char *fname, ChannelInfo *chans );

  /* ------------- WRITER --------------------- */
  void write_double_dblpp_ascii(FILE *out, const double **d, int xdim, int ydim, OptArgList *opts);
  void write_double_dblpp_ascii_file(const char *fname, const double **d, int xdim, int ydim, OptArgList *opts);

  void write_double_dblp_ascii(FILE *out, const double *v, int n);
  void write_double_dblp_ascii_file(const char *fname, const double *v, int n);

  void write_int_dblp_ascii(FILE *out, const int *v, int n);
  void write_int_dblp_ascii_file(const char *fname, const int *v, int n);


#ifdef __cplusplus
}
#endif

#endif /* IO_H */
