/**\file helper.h
 * \brief Helper functions.
 *
 * \defgroup helper Helper functions
 *\{
 * \defgroup helpermath Helper Math functions
 * \defgroup helperorg Organizing Helper functions
 * \defgroup helperplot Plotting Functions
 *\}
 */
#ifndef HELPER_H
# define HELPER_H

/* #define HAVE_PLOTLIB */
/*#define HAVE_MATLAB*/

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
//#define DEBUG
#ifdef DEBUG
/*#define dprintf(...) fprintf(stderr, ## __VA_ARGS__)*/
#define dprintf(...) do{																\
		fprintf(stderr, ESCAPE_FGYELLOW_STR ESCAPE_BRIGHT_STR "%s (%i), %s(): ", __FILE__, __LINE__, __FUNCTION__); \
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

#define oprintf(...) do{ fprintf(stderr, ESCAPE_FGCYAN_STR ESCAPE_BRIGHT_STR \
											"%s (%i), %s() --> ",						\
											__FILE__, __LINE__, __FUNCTION__);		\
	 fprintf(stderr, ## __VA_ARGS__);												\
	 fprintf(stderr, ESCAPE_RESET_STR); } while(0)

#ifdef __cplusplus
extern "C" {
#endif

  /* ---------------------------------------------------------------------------- 
	  -- Helper functions                                                       -- 
	  ---------------------------------------------------------------------------- */
  /**\addtogroup helperorg
	*\ingroup helper
	*\{*/

  int     v_printf(int v, char *format, ...);
  int     vprint_vector(const char* name, double *v, int n);
  void    errormsg(int err_no, int fatal);
  size_t  ffread(void *ptr, size_t size, size_t nmemb, FILE *stream);

  void    swap_bytes(void *ptr, int nmemb);
  void    wswap(void *ptr, int nmemb, int flag);
  int     is_little_endian();

  void    qsort_int_index( int *idx_idx, const int *idx, int n );
  int     compare_ints (const int *a, const int *b);

  void    copy_modeldata(const ModelData *m1, ModelData *m2);


  /* constructors */
  EEGdata* init_eegdata(int nbchan, int nsamples, int nmarkers);
  EEGdata_trials* init_eegdata_trials(int nbtrials, int markers_per_trial, int nbchan, int nbsamples);
  WarpPath* init_warppath(int J, int K);

  void      reset_warppath(WarpPath *P, int J, int K);

  /* destructors */
  void    free_modeldata(ModelData *m);
  void    free_warppath(WarpPath *p);
  void    free_eegdata(EEGdata *eeg);
  void    free_eegdata_trials(EEGdata_trials *eeg);

  /* printing-functions */
  void    print_modeldata(FILE *out, const ModelData *m);
  void    print_denoisingparameters(FILE *out, const DenoisingParameters *p);
  void    print_timewarpparameters(FILE *out, const TimewarpParameters *p);
  void    print_eegdata_trials(FILE *out, const EEGdata_trials *eeg);

  double** copy_double_ptrptr(const double **s, int N, int n);


/**\}*/


/**\addtogroup helperplot
 *\ingroup helper
 *\{*/
#ifdef HAVE_MATLAB
Engine*  ml_init(void);
int  ml_close(Engine *m);
void ml_plot(Engine *matlab, const double *r, const double *v, int n, const char* color, int new);
void ml_plot_path(Engine *matlab, const int *path, int K);
void ml_wait(Engine *matlab);
#endif

/**\}*/

/* ---------------------------------------------------------------------------- 
   -- Gobals                                                                 -- 
   ---------------------------------------------------------------------------- */
int *verbosity; 
#ifdef HAVE_MATLAB
Engine *matlab;
#endif

#ifdef __cplusplus
}
#endif

#endif
