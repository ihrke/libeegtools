/**\file denoising.h
 * \brief Contains functions used for signal denoising.
 *
 * Especially wavelet-based denoising. 
 */
/** Groups.
 * \defgroup denoising Denoising Routines 
 * \{
 *    \defgroup wavelet Wavelet-Denoising
 *    \{
 *       \defgroup thresholding Thresholding Functions
 *       \defgroup select_thresh Threshold selection Functions
 *       \defgroup sigext Signal Extension Functions
 *    \}
 *    \defgroup robust_filtering Robust-Filtering
 * \}
 *
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
   -- Robust filtering methods  --
   ------------------------------ */

/** \addtogroup robust_filtering
 * \ingroup denoising
 *\{
 */
void eeg_filter_running_median(EEGdata *s, int win);
double* running_median(double *d, int n, int win);
void eeg_filter_weighted_running_median(EEGdata *s, int win);
double* weighted_running_median(double *d, int n, int win, 
				double(*dist)(double,double));
double weighted_median_from_unsorted(const double *d, const double *w, int n);
double dist_euclidean(double x, double y);
/** \} */


/* ---------------------------------------------------------------------------- 
   -- Wavelet-based Denoising routines                                                     -- 
   ---------------------------------------------------------------------------- */

/** \addtogroup thresholding
 * \ingroup wavelet
 *\{
 */
double eta_s(double d, double lambda);
double eta_h(double d, double lambda);
/** \} */

/** \addtogroup select_thresh
 *\ingroup wavelet
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
   -- Signal extension routines                                              -- 
   ---------------------------------------------------------------------------- */

/** \addtogroup sigext
 * \ingroup wavelet
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
 *
 * \{
 */
double* sigext_zeros(double *data, int ns, int n);
double* sigext_zerosr(double *data, int ns, int n);
double* sigext_sym(double *data, int ns, int n);
double* sigext_smooth(double *data, int ns, int n);
/** \} */

#endif /* - DENOISING_H - */
