/* **************************************************************************
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

/**\file helper.h
 * \brief Helper functions.
 */
#ifndef HELPER_H
# define HELPER_H

#include <stdarg.h>
#include <string.h> /* memcpy */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <limits.h>
#include <float.h>
#include <stdint.h>

#define _WITH_ANSI_ESCAPE_CODES /* necessary for the esc-codes to show up */
#include "escape_codes.h"

#include "definitions.h"

#ifdef HAVE_MATLAB
#ifndef MATLAB_STARTUP_CMD
#define MATLAB_STARTUP_CMD "/usr/local/matlab7/bin/matlab -nosplash"   
/* #define MATLAB_STARTUP_CMD "/usr/nld/matlab-17/bin/matlab -nosplash" */
#endif
#define ML(code) code /* this macro deletes all ML() statements if
								 HAVE_MATLAB is not defined */
#include "engine.h"
#else
#define ML(code)
#endif

/** custom error stuff */
#define ERR_GSL  -1
#define ERR_IO   -2
#define ERR_PLOT -3
#define ERR_MEM  -4
#define ERR_MATLAB -5
#define ERR_ENDIAN -6
#define ERR_PARSEMAT -7
#define ERR_FATAL 1
#define ERR_NOFATAL 0


/** debugging macros */
#ifdef DEBUG
/*#define dprintf(...) fprintf(stderr, ## __VA_ARGS__)*/
#define dprintf(...) do{																\
		fprintf(stderr, ESCAPE_FGYELLOW_STR ESCAPE_BRIGHT_STR "%s (%i), %s(): ", \
				  __FILE__, __LINE__, __FUNCTION__);								\
		fprintf(stderr, ## __VA_ARGS__);												\
		fprintf(stderr, ESCAPE_RESET_STR);											\
	 } while(0)

#define massert(x, text, ...) \
  do{																							\
	 if(x){																					\
		fprintf(stderr, ESCAPE_FGRED_STR "Assertion failed: " #text, ## __VA_ARGS__); \
		fprintf(stderr,  ESCAPE_RESET_STR);											\
	 }																							\
  } while(0)
#else
#define dprintf(...)
#define massert(...)	    
#endif

#define errprintf(...) do{ fprintf(stderr, ESCAPE_FGRED_STR ESCAPE_BOLD_STR\
											  ESCAPE_BGYELLOW_STR "ERROR: %s (%i), %s(): ", \
											  __FILE__, __LINE__, __FUNCTION__);	\
	 fprintf(stderr, ## __VA_ARGS__);												\
	 fprintf(stderr, ESCAPE_RESET_STR); } while(0)

#define warnprintf(...) do{ fprintf(stderr, ESCAPE_FGWHITE_STR ESCAPE_BOLD_STR\
											  ESCAPE_BGBLUE_STR "WARNING: %s (%i), %s(): ", \
											  __FILE__, __LINE__, __FUNCTION__);	\
	 fprintf(stderr, ## __VA_ARGS__);												\
	 fprintf(stderr, ESCAPE_RESET_STR); } while(0)

#define oprintf(...) do{ fprintf(stderr, ESCAPE_FGCYAN_STR ESCAPE_BRIGHT_STR \
											"%s (%i), %s() --> ",						\
											__FILE__, __LINE__, __FUNCTION__);		\
	 fprintf(stderr, ## __VA_ARGS__);												\
	 fprintf(stderr, ESCAPE_RESET_STR); } while(0)

#ifdef __cplusplus
extern "C" {
#endif

  /* ---------------------------------------------------------------------------- 
	  -- Helper functions (org)                                                       -- 
	  ---------------------------------------------------------------------------- */
  /**\addtogroup helperorg
	*\{*/
  void     errormsg(int err_no, int fatal);

  void     swap_bytes(void *ptr, int nmemb);
  void     wswap(void *ptr, int nmemb, int flag);
  int      is_little_endian();

  void     qsort_int_index( int *idx_idx, const int *idx, int n );
  int      compare_ints (const void *a, const void *b);

  int      randint( int from, int to );  

  int      strcount( const char *s, char c );

  double** copy_double_ptrptr(const double **s, int N, int n);
  /**\}*/

  /**\addtogroup helperio
	*\{*/
  /* ---------------------------------------------------------------------------- 
	  -- Helper functions (IO)                                                  -- 
	  ---------------------------------------------------------------------------- */
  int     v_printf(int v, char *format, ...);
  int     vprint_vector(const char* name, double *v, int n); 

  size_t  ffread(void *ptr, size_t size, size_t nmemb, FILE *stream);
  size_t  ffwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream);

  int     stream_count_char( FILE* f, char c );

  /**\}*/

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

  /* ---------------------------------------------------------------------------- 
	  -- Gobals                                                                 -- 
	  ---------------------------------------------------------------------------- */

  /** \addtogroup progressbar 
		\{  
  */
  static ProgressBarStatus progress_status;
  void   progressbar_rotating( int flag, int num );
  /** \} */

  static int *verbosity; 

#ifdef __cplusplus
}
#endif

#endif
