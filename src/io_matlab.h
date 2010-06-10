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

/**\file io_matlab.h
 \brief \ref status_inprogress io_matlab MATLAB-input-output for (EEG)-data.

 eeglab/matlab file reader functions if matio library is installed 
	
 */
#ifndef IO_MATLAB_H
# define IO_MATLAB_H

#ifdef __cplusplus
extern "C" {
#endif

#include "array.h"
#include "eeg.h"
  
#ifdef MATIO

  /* -------------- READER ---------------- */
  EEG* read_eeglab_file( const char *file );
  Array* read_array_matlab( const char *file, const char *varname ); 

  /* -------------- WRITER ---------------- */
  int write_eeglab_file( EEG* eeg, const char *file );
  int write_array_matlab( const Array *a, const char *varname, 
								  const char *file, bool append );
#endif



#ifdef __cplusplus
}
#endif

#endif /* IO_MATLAB_H */
