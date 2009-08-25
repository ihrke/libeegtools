%{
#include "definitions.h"
#include "helper.h"
%}


typedef struct{
  unsigned int nbchan;    /**< number of channels */
  unsigned int n;         /**< number of samples */
  double **d;             /**< data */
  unsigned long *markers; /**< stimulus onset, response, etc
									  in sampling points; nmarkers long */
  unsigned int nmarkers;
} EEGdata;

typedef struct{
  EEGdata **data;          /**< EEG-data for all  trials (n EEGdata-struct-ptr) */
  unsigned int ntrials;    /**< number of trials = dim(eeg) */  
  unsigned int nmarkers_per_trial;
  unsigned long **markers; /**< stimulus onset, response, etc
										in sampling points; dims are ntrials x nmarkers */
  unsigned int nsamples;   /**< number of samples in each data */
  double *times;           /**< times array; data->n long */
  double sampling_rate;    /**< in Hz */
} EEGdata_trials;


EEGdata*        init_eegdata(int nbchan, int nsamples, int nmarkers);
void            reset_eegdata( EEGdata* eeg );
EEGdata_trials* init_eegdata_trials(int nbtrials, int markers_per_trial, 
												int nbchan, int nbsamples, double *times);


/* destructors */
void    free_eegdata(EEGdata *eeg);
void    free_eegdata_trials(EEGdata_trials *eeg);

/* convenience functions for structs */
int       eegdata_cmp_settings( const EEGdata *s1, const EEGdata *s2 );

EEGdata_trials* clone_eegdata_trials( const EEGdata_trials *source );
EEGdata*        clone_eegdata( const EEGdata *source );



%extend EEGdata_trials{
  EEGdata_trials(int nbtrials, int markers_per_trial, 
					  int nbchan, int nbsamples ){
	 return init_eegdata_trials( nbtrials,  markers_per_trial, 
										  nbchan, nbsamples, NULL);
  }

  EEGdata_trials* clone(){
	 return clone_eegdata_trials( self );
  }

  EEGdata* get_trial( int i ){
	 if( i>=0 || i<self->ntrials ){
		return self->data[i];
	 } else {
		return NULL;
	 }
  }

  PyObject* get_times(){
	 PyObject *l;
	 l = doubleptr_to_pythonlist( self->times, self->nsamples );
	 return l;
  }

  char* __str__(){
	 static char buf[512];
	 char *tmp;
	 
	 sprintf(buf, "EEGdata_trials:\n"  
				"      ntrials=%i\n"
				"      nmarkers_per_trial=%i\n", 
				self->ntrials, self->nmarkers_per_trial);
	 tmp=&(buf[strlen(buf)]);
	 if( self->nmarkers_per_trial>0 ){
		sprintf(tmp, "        markers[0][0]=%ld\n"
				  "        markers[0][nm-1]=%ld\n", 
				  self->markers[0][0],self->markers[0][self->nmarkers_per_trial-1]);
	 }
	 tmp=&(buf[strlen(buf)]);
	 sprintf(tmp, "      sampling_rate=%f\n", self->sampling_rate);
	 tmp=&(buf[strlen(buf)]);
	 if( self->times!=NULL ){
		sprintf(tmp, "      times[0]=%f\n"
				  "      times[n-1]=%f\n", 
				  self->times[0], self->times[self->nsamples-1]);
	 } else {
		sprintf(tmp, "      times[0]=<NULL>\n");
	 }
	 tmp=&(buf[strlen(buf)]);
	 if( self->data[0] && self->data[0]->d[0] && self->data[0]->d[0][0] ){
		sprintf(tmp, "      data[0]->d[0][0]=%f\n"
				  "      data[0]->d[0][n-1]=%f\n", 
				  self->data[0]->d[0][0], self->data[0]->d[0][self->nsamples-1] );
	 } else {
		sprintf(tmp, "      data[0]->d[0][0]=<NULL>\n" );
	 }
	 return buf;
  }
};

%extend EEGdata {
  EEGdata( int nbchan, int nsamples, int nmarkers ){
	 return init_eegdata( nbchan, nsamples, nmarkers);
  }
  void reset(){
	 reset_eegdata(self);
  }
 
  EEGdata* clone(){
	 return clone_eegdata( self );
  }

  PyObject* get_chan( int chan ){
	 if( chan<0 || chan>=self->nbchan ){
		PyErr_SetString(PyExc_TypeError,"channel not available");
		return NULL;
	 }
	 PyObject *l = PyList_New( self->n );
	 int i;
	 for( i=0; i<self->n; i++ ){
		PyList_SetItem( l, i, PyFloat_FromDouble( self->d[chan][i] ));
	 }
	 return l;
  }

  PyObject* get_data( ){
	 PyObject *o, *l = PyList_New( self->nbchan );
	 int i,j;
	 for( i=0; i<self->nbchan; i++ ){
		o = PyList_New( self->n );	  
		for( j=0; j<self->n; j++ ){
		  PyList_SetItem( o, j, PyFloat_FromDouble( self->d[i][j] ));
		}
		PyList_SetItem( l, i, o );
	 }
	 return l;
  }

  char* __str__(){
	 static char buf[256];
	 char *tmp;
	 
	 sprintf(buf, "EEGdata:\n"
				"      nbchan  =%i\n"
				"      n       =%i\n"
				"      nmarkers=%i\n", self->nbchan, self->n, self->nmarkers);
	 tmp=&(buf[strlen(buf)]);
	 if( self->nmarkers>0 ){
		sprintf(tmp, 
				  "        markers[0   ]=%ld\n"
				  "        markers[nm-1]=%ld\n", 
				  self->markers[0],self->markers[self->nmarkers-1]);
	 }
	 if( self->n>0 ){
		sprintf(tmp, 
				  "      d[0][0  ]=%f\n"
				  "      d[0][n-1]=%f\n", 
				  self->d[0][0],self->d[0][self->n-1]);
	 }
	 return buf;
  }
};
