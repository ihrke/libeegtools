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
 * \brief \ref status_stable Helper functions.
 */
#ifndef HELPER_H
# define HELPER_H

#include "definitions.h"
#include "escape_codes.h"
#include <stdio.h>


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

#define NOSTRING " \b" 

#define PROGRESSBAR_NUMCOLS        80 /**< number of columns for progress bar */
#define PROGRESSBAR_INIT           0  /**< initialize flag for progress bar   */
#define PROGRESSBAR_CONTINUE_SHORT 1  /**< continue flag for progress bar, minor step */
#define PROGRESSBAR_CONTINUE_LONG  2  /**< continue flag for progress bar, major step */
#define PROGRESSBAR_FINISH         3  /**< cleanup */

/** \brief swap two variables.
 */
#define SWAPT( Type, a, b ) { Type tmp; tmp=(a); (a)=(b); (b)=(tmp); }


/*-----------------------------------------------------------
  - Progress Bar
  ---------------------------------------------------------*/

/** \brief status of progress bar.
	 This struct is used as an external variable to
	 keep track of the progress bar. You only need
	 to be concerned about this, if you want to 
	 write your own progressbar-function.
 */
typedef struct{
  ulonglong max_progress;
  ulonglong cur_progress;
  ulonglong prev_progress;
} ProgressBarStatus;

/** \brief A progressbar-function.
 */ 
typedef void    (*ProgressBarFunction)          (int,ulonglong);

  /** \brief convenience memory alloc macro. 

		use as follows:
		\code
		double *a;
		MALLOC( a, 100, double );
		double **b;
		MALLOC( b, 100, double* ); 
		\endcode
  */
#define MALLOC( ptr, nmemb, type )						\
  ptr = (type*)malloc(nmemb*sizeof(type));			\
  if( ptr==NULL ){											\
	 errprintf("Could not allocate memory\n");		\
  }
  
/** debugging macros */
#ifdef DEBUG
/*#define dprintf(...) fprintf(stderr, ## __VA_ARGS__)*/

/** \brief debug-print.

	 This is only printed if the library is compiled with
	 the DEBUG flag.
*/
#define dprintf(...) do{																\
		fprintf(stderr, ESCAPE_FGYELLOW_STR ESCAPE_BRIGHT_STR "%s (%i), %s(): ", \
				  __FILE__, __LINE__, __FUNCTION__);								\
		fprintf(stderr, ## __VA_ARGS__);												\
		fprintf(stderr, NOSTRING ESCAPE_RESET_STR);											\
	 } while(0)

#define massert(x, text, ...) \
  do{																							\
	 if(x){																					\
		fprintf(stderr, ESCAPE_FGRED_STR "Assertion failed: " #text, ## __VA_ARGS__); \
		fprintf(stderr, NOSTRING ESCAPE_RESET_STR);											\
	 }																							\
  } while(0)
#else
#define dprintf(...)
#define massert(...)	    
#endif

/** \brief print error.

	 Report an Error.
*/
#define errprintf(...) do{ fprintf(stderr, ESCAPE_FGRED_STR ESCAPE_BOLD_STR\
											  ESCAPE_BGYELLOW_STR "ERROR: %s (%i), %s(): ", \
											  __FILE__, __LINE__, __FUNCTION__);	\
	 fprintf(stderr, ## __VA_ARGS__);												\
	 fprintf(stderr, NOSTRING ESCAPE_RESET_STR); } while(0)
/** \brief print warning.
*/
#define warnprintf(...) do{ fprintf(stderr, ESCAPE_FGWHITE_STR ESCAPE_BOLD_STR \
												ESCAPE_BGBLUE_STR "WARNING: %s (%i), %s(): ", \
												__FILE__, __LINE__, __FUNCTION__);	\
	 fprintf(stderr, ## __VA_ARGS__);												\
	 fprintf(stderr, NOSTRING ESCAPE_RESET_STR); } while(0)

#define oprintf(...) do{ fprintf(stderr, ESCAPE_FGCYAN_STR ESCAPE_BRIGHT_STR \
											"%s (%i), %s() --> ",						\
											__FILE__, __LINE__, __FUNCTION__);		\
	 fprintf(stderr, ## __VA_ARGS__);												\
	 fprintf(stderr, NOSTRING ESCAPE_RESET_STR); } while(0)

#ifdef __cplusplus
extern "C" {
#endif

  /* ---------------------------------------------------------------------------- 
	  -- Helper functions (org)                                                       -- 
	  ---------------------------------------------------------------------------- */
  void     errormsg(int err_no, int fatal);

  void     swap_bytes(void *ptr, int nmemb);
  void     wswap(void *ptr, int nmemb, int flag);
  int      is_little_endian();

  bool  isin_intarray( const int *a, int n, int val );

  void     qsort_int_index( int *idx_idx, const int *idx, int n );
  int      compare_ints (const void *a, const void *b);

  int      randint( int from, int to );  

  int      strcount( const char *s, char c );

  int      safer_free( void *p );

  char*    create_string( const char *string );

  double** copy_double_ptrptr(const double **s, int N, int n);
  void     string_strip_blanks( char *s );

  /* ---------------------------------------------------------------------------- 
	  -- Helper functions (IO)                                                  -- 
	  ---------------------------------------------------------------------------- */
  int     vprint_vector(const char* name, double *v, int n); 

  size_t  ffread(void *ptr, size_t size, size_t nmemb, FILE *stream);
  size_t  ffwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream);

  int     stream_count_char( FILE* f, char c );




  /* ---------------------------------------------------------------------------- 
	  -- Gobals                                                                 -- 
	  ---------------------------------------------------------------------------- */

  extern ProgressBarStatus progress_status;
  void   progressbar_rotating( int flag, ulonglong num );

#ifdef __cplusplus
}
#endif

#endif
