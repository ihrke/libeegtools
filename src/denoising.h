/** \mainpage 
 *
 * These files contain signal denoising and timewarping functions.
 *
 * Synopsis:
 * - ml_*.c -- contain MATLAB-wrapper MEX-files that can be called
 *             directly from MATLAB
 *    - ml_timewarp.c 
 *    - ml_denoise.c
 *    - ml_warpavg.c
 * - t_*.c -- contain test programs for the parts of the algorithm.
 *    - t_denoise.c
 *    - t_denoisestat.c
 *    - t_timewarpstat.c
 *    - t_warpavg.c
 * - helper.h, helper.c contain helper math, print, plot functions
 * etc.
 * - denoising.h, denoising.c - core functionality
 * - argspec.h argspec.c - cmd-lind parsing for the t_* programs
 * - directory \link test test \endlink contains checklib unit tests
 *
 * \page history History.
 *\code
 * {{HISTORY}}
 *\endcode
 */

/**\file denoising.h
 * \brief Contains functions used for signal denoising.
 *
 * Especially wavelet-based denoising. 
 */
#ifndef DENOISING_H
# define DENOISING_H

#include <gsl/gsl_wavelet.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdarg.h>
#include <string.h> /* memcpy */
#include "helper.h"
#include "mathadd.h"


/* ------------------------------ 
   -- Other filtering methods  --
   ------------------------------ */

void eeg_filter_running_median(EEGdata *s, int win);
double* running_median(double *d, int n, int win);
void eeg_filter_weighted_running_median(EEGdata *s, int win);
double* weighted_running_median(double *d, int n, int win, 
				double(*dist)(double,double));
double weighted_median_from_unsorted(const double *d, const double *w, int n);
double dist_euclidean(double x, double y);

/* ---------------------------------------------------------------------------- 
   -- Denoising routines                                                     -- 
   ---------------------------------------------------------------------------- */
/** Groups.
 * \defgroup denoising Denoising Routines 
 * \{
 * \defgroup thresholding Thresholding Functions
 * \defgroup select_thresh Threshold selection Functions
 * \defgroup sigext Signal Extension Functions.
 * \}
 * \defgroup timewarp Timewarping Functions
 *
 */

/** Implements soft-thresholding.
 * Formula: \f$ \eta_s(\lambda, w)=\f$
 * \ingroup thresholding
 */
double eta_s(double d, double lambda);
/** Implements soft-thresholding.
 * Formula: \f$ \eta_s(\lambda, w)=\f$
 * \ingroup thresholding
 */
double eta_h(double d, double lambda);


/** \addtogroup select_thresh
 *\{
 */
int generic_denoising(double *data, int n, int L, 
		      double(*threshfct)(const double*, int), 
		      double(*etafct)(double,double));
int extend_and_denoise(double *data, int n, int L, 
		       double(*threshfct)(const double*, int), 
		       double(*etafct)(double,double), 
		       double*(*sigextfct)(double*, int, int));
double translation_invariant_thresholding(const double *data, int n);
double conventional_thresholding(const double *data, int n);
double sureshrink(const double *data, int n);
double heuristic_sure(const double *data, int n);
/** \} */

/* ---------------------------------------------------------------------------- 
   -- Merit Measures                                                         -- 
   ---------------------------------------------------------------------------- */
double rmse(const double *r, const double *d, int n);
double snr (const double *r, const double *d, int n);


/* ---------------------------------------------------------------------------- 
   -- Signal extension routines                                              -- 
   ---------------------------------------------------------------------------- */

/** \addtogroup sigext
 * Signal extension schemes to extend signal of length n to length 2^j
 *    with 2^j being the closest power of 2 to n.
 * The extension functions return a pointer to the former data[0],
 * because this is where the unextended signal began;\n
 * Example:
 * \code
 * sigext([1 2 3 - - - -]) -> [0 0 1 2 3 0 0] 
 *                                 ^ ptr
 * \endcode
 *  Assumptions (not for full generality!):
 * -# ns <= n
 * -# n <= 2*ns
 * \ingroup denoising
 * \{
 */
double* sigext_zeros(double *data, int ns, int n);
double* sigext_zerosr(double *data, int ns, int n);
double* sigext_sym(double *data, int ns, int n);
double* sigext_smooth(double *data, int ns, int n);
/** \} */

#endif /* - DENOISING_H - */
