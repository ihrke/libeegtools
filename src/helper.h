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
#include "definitions.h"

#ifdef HAVE_PLOTLIB
#include <plot.h>
#endif

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

#define ERR_GSL  -1
#define ERR_IO   -2
#define ERR_PLOT -3
#define ERR_MEM  -4
#define ERR_MATLAB -5
#define SCREENX  1280
#define SCREENY  800
#define YSCALE   120
#define DEBUG

#ifdef DEBUG
#define dprintf(...) fprintf(stderr, ## __VA_ARGS__)
#define massert(x, text, ...) \
if(x) fprintf(stderr, "Assertion failed: " #text, ## __VA_ARGS__)
#else
#define dprintf(...)
#define massert(...)	    
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

void    copy_modeldata(const ModelData *m1, ModelData *m2);
void    free_modeldata(ModelData *m);
void    free_warppath(WarpPath *p);
void    print_modeldata(FILE *out, const ModelData *m);
void    print_denoisingparameters(FILE *out, const DenoisingParameters *p);
void    print_timewarpparameters(FILE *out, const TimewarpParameters *p);
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

#ifdef HAVE_PLOTLIB
int    plot_init(double xmin, double xmax, double ymin, double ymax, char *plotdev);
void   plot_trace(const double *times, const double *d, int n, const char *color);
void   plot_coordsys(double xmin, double xmax, double ymin, double ymax, double *labels, int numlab);
void   plot_close(int handle);
#endif
/**\}*/

/* ---------------------------------------------------------------------------- 
   -- Gobals                                                                 -- 
   ---------------------------------------------------------------------------- */
int *verbosity; 
#ifdef HAVE_MATLAB
Engine *matlab;
#endif

#endif
