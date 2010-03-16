%{
#include "filter.h"
#include "eeg.h"
#include "wavelet.h"
#include "definitions.h"
%}


typedef double  (*ThresholdSelectionFunction)   (const double*,int);
typedef double  (*ThresholdFunction)            (double,double);
typedef double* (*SignalExtensionFunction)      (double*,int,int);

%apply (double* INPLACE_ARRAY1, int DIM1) {(double* d, int n)}
double* running_median(double *d, int n, int win);
double* moving_average(double *d, int n, int win);
double* weighted_running_median( double *d, int n, int win, 
											PointDistanceFunction dist);


EEG*    eeg_filter_running_median( EEG *eeg, int win, bool alloc );
EEG*    eeg_filter_weighted_running_median( EEG *eeg, int win, bool alloc );

EEG* eeg_filter_fidlib( EEG *eeg, const char *spec, bool alloc );


%typemap(in) (double *data, int n){
  $1 = python_list_to_doubleptr($input);
  $2 = PyList_Size($input);
}

%callback("%s_cb");
double translation_invariant_thresholding( const double *data, int n );
double conventional_thresholding         ( const double *data, int n );
double sureshrink                        ( const double *data, int n );
double heuristic_sure                    ( const double *data, int n );




double* sigext_zeros(double *data, int ns, int n);
double* sigext_zerosr(double *data, int ns, int n);
double* sigext_sym(double *data, int ns, int n);
double* sigext_smooth(double *data, int ns, int n);
%nocallback;


/* %pythoncode */
/* %{ */
/* def wavelet_denoising(eeg, L=4, thresholdselection=heuristic_sure_cb, */
/*                       thresholdfct=eta_s_cb, signalextension=sigext_sym_cb, */
/*                       trials=None): */
/*     if type(eeg).__name__=='EEGdata': */
/*         print "This is an EEGdata struct"; */
/*         neweeg = eeg.clone(); */
/*         eeg_wavelet_denoise( neweeg, L, thresholdselection, thresholdfct, signalextension ); */
/*         return neweeg */
/*     elif type(eeg).__name__=="EEGdata_trials": */
/*         print "This is an EEGdata_trials struct" */
/*         if not trials: */
/*             trials = range(0,eeg.ntrials); */
/*         neweeg = eeg.clone(); */
/*         for i in trials: */
/*             x=neweeg.get_trial(i); */
/*             eeg_wavelet_denoise( x, L, thresholdselection, thresholdfct, signalextension ); */
/*         return neweeg */
/*     else: */
/*         print "Plotting does not support '%s'"%(type(eeg).__name__) */

/* %} */
