%{
#include "definitions.h"
#include "eeg.h"
#include "reader.h"
%}

#define EEG_CLONE_ALL        0     /**< clone everything */
#define EEG_CLONE_NODATA     2<<1  /**< clone everything except 
												  the eeg->data field */
#define EEG_CLONE_NOMARKERS  2<<2  /**< clone everything except
												  the eeg->markers, eeg->nmarkers
												  eeg->marker_labels fields */
#define EEG_CLONE_NOCHANINFO 2<<3  /**< clone everything except
												  the eeg->chaninfo field */


  /*-----------------------------------------------------------
	 - EEG-CHANNELS -
	 ---------------------------------------------------------*/

  typedef struct{
	 int    num;
	 int    num_chans;
	 char   label[MAX_LABEL_LENGTH];
	 double x;
	 double y;
	 double z;
  } ChannelInfo;


  /*-----------------------------------------------------------
	 - EEG-Data -
	 ---------------------------------------------------------*/

  typedef struct{
	 char         *filename;
	 char         *comment;
	 unsigned int nbchan;  /**< number of channels */
	 unsigned int ntrials; /**< number of trials = dim(eeg) */ 
	 unsigned int n;       /**< number of samples */
    double       sampling_rate; /**< in Hz */
	 double       *times;  /**< times array; n-long */
	 ChannelInfo  *chaninfo; /**< location and other information
										 about the channels */

	 double       ***data; /**< channels x trials x samples */

	 unsigned int *nmarkers; /**< number of markers for each trial */
	 unsigned int **markers; /**< trials x nmarkers[i] */
	 char       ***marker_labels; /**< trials x nmarkers[i] x length(label) gives  a 
											  label for each marker */
  } EEG;

EEG* read_eeglab_file( const char *file );  
ChannelInfo* read_chaninfo_ced( const char *fname, ChannelInfo *chans );

%extend ChannelInfo{
  ChannelInfo __getitem__(int i){
	 return self[i];
  }
}

%extend EEG{
  EEG(int nbchan, int ntrials, int nsamples ){
	 return eeg_init( nbchan, ntrials, nsamples );
  }

  EEG* clone( int clone_flags ){
	 return eeg_clone( self, clone_flags );
  }

  PyObject* get_data(){
	 npy_intp dims[3];
	 dims[0] = self->nbchan;
	 dims[1] = self->ntrials;
	 dims[2] = self->n;
	 PyArrayObject* a;
	 a=(PyArrayObject *)PyArray_SimpleNewFromData( 3, dims, NPY_DOUBLE, self->data );
	 return (PyObject*)a;
  }

  PyObject* get_times(){
	 PyArrayObject *l;
	 npy_intp dims=self->n;
	 if( self->times == NULL ) {
		Py_INCREF(Py_None);
		return Py_None;
	 }
	 
	 l = (PyArrayObject *)PyArray_SimpleNewFromData( 1, &dims, 
																	 NPY_DOUBLE, self->times );
	 return (PyObject*)l;
  }

  char* __str__(){
	 static char *buf;
	 buf = eeg_sprint(NULL, self, 2 );
	 return buf;
  }
};

