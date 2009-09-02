%{
#include <libeegtools/denoising.h>
%}

%typemap(in) (double *d, int n){
  $1 = python_list_to_doubleptr($input);
  $2 = PyList_Size($input);
 }
%typemap(argout) (double*){
  int i; 
  int l=PyList_Size($input);
  PyObject *o;
  $result = PyList_New( 0 );
  for( i=0; i<l; i++ ){
	 o = PyFloat_FromDouble( $1[i] );
	 PyList_Append( $result, o );
  }
 }

typedef double  (*ThresholdSelectionFunction)   (const double*,int);
typedef double  (*ThresholdFunction)            (double,double);
typedef double* (*SignalExtensionFunction)      (double*,int,int);


double* running_median(double *d, int n, int win);
double* moving_average(double *s, int n, int win);


void eeg_filter_running_median(EEGdata *eeg, int win);
void eeg_filter_fidlib( EEGdata *eeg, double sampling_rate, const char *spec ); 
void eeg_wavelet_denoise( EEGdata *eeg, int L, 
								  ThresholdSelectionFunction threshfct,
								  ThresholdFunction etafct, 
								  SignalExtensionFunction sigextfct );

%typemap(in) (double *data, int n){
  $1 = python_list_to_doubleptr($input);
  $2 = PyList_Size($input);
}
int extend_and_denoise  ( double *data, int n, int L, 
								  ThresholdSelectionFunction threshfct,
								  ThresholdFunction etafct, 
								  SignalExtensionFunction sigextfct );

%callback("%s_cb");
double translation_invariant_thresholding( const double *data, int n );
double conventional_thresholding         ( const double *data, int n );
double sureshrink                        ( const double *data, int n );
double heuristic_sure                    ( const double *data, int n );

double eta_s(double d, double lambda);
double eta_h(double d, double lambda);

double* sigext_zeros(double *data, int ns, int n);
double* sigext_zerosr(double *data, int ns, int n);
double* sigext_sym(double *data, int ns, int n);
double* sigext_smooth(double *data, int ns, int n);
%nocallback;


%pythoncode
%{
def wavelet_denoising(eeg, L=4, thresholdselection=heuristic_sure_cb,
                      thresholdfct=eta_s_cb, signalextension=sigext_sym_cb,
                      trials=None):
    if type(eeg).__name__=='EEGdata':
        print "This is an EEGdata struct";
        neweeg = eeg.clone();
        eeg_wavelet_denoise( neweeg, L, thresholdselection, thresholdfct, signalextension );
        return neweeg
    elif type(eeg).__name__=="EEGdata_trials":
        print "This is an EEGdata_trials struct"
        if not trials:
            trials = range(0,eeg.ntrials);
        neweeg = eeg.clone();
        for i in trials:
            x=neweeg.get_trial(i);
            eeg_wavelet_denoise( x, L, thresholdselection, thresholdfct, signalextension );
        return neweeg
    else:
        print "Plotting does not support '%s'"%(type(eeg).__name__)

%}
