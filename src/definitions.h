/** \mainpage 
 *
 * This is libEEGTools. Please explore the functionality via the 
 * Modules page (link above). 
 *
 * The official Homepage is http://pinguin.uni-psych.gwdg.de/~ihrke/wiki/index.php/LibEEGTools 
 *
 * \section synopsis Synopsis:
 * - denoising.h, denoising.c - core functionality
 * - averaging.h, averaging.c - averaging function
 * - clustering.h clustering.c - clustering functions
 * - reader.c, reader.h - get data
 * - helper.h, helper.c contain helper math, print, plot functions etc.
 * - mathadd.h, mathadd.c - mathematical helper
 *
 * \section notes Notes
 * - Functions starting with eeg_*() take an EEGdata - struct as argument 
 *   and the algorithm is carried out independantly for each channel within
 *   the struct.
 * - Functions that start with eegtrials_*() take an EEGdata_trials - struct
 *   as main argument and the algorithm is carried out independantly on the
 *   trial segments (usually the corresponding eeg_*() function is called for
 *   each trial-segment)
 *
 *
 *
 * History:
\code
 * $Log$
 * Revision 1.4  2008/05/29 12:40:37  mihrke
 * version 0.2
 *
 * Revision 1.3  2008/03/14 15:44:27  mihrke
 * distribution 0.15
 *
\endcode
 */

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

typedef struct{
  unsigned int nbchan; /** number of channels */
  unsigned int n; /** number of samples */
  double **d; /** data */
  unsigned long *markers; /** stimulus onset, response, etc
										in sampling points; nmarkers long */
  unsigned int nmarkers;
} EEGdata;

typedef struct{
  EEGdata **data; /** EEG-data for all  trials (n EEGdata-struct-ptr) */
  unsigned int ntrials; /** number of trials = dim(eeg) */  
  unsigned int nmarkers_per_trial;
  unsigned long **markers; /** stimulus onset, response, etc
										 in sampling points; dims are ntrials x nmarkers */
  double *times; /** times array; data->n long */
} EEGdata_trials;

/**\ingroup clustering */
typedef struct{
  int **clust; /** indices for the trials in the cluster (Kxn)*/
  int K;  /** number of clusters */
  int *n; /** number of trials in each of the K clusters */
} Clusters;


/** \ingroup denoising */
typedef struct{
	  int L;
	  double (*cleanfct)(const double*, int); /** cleaning function to use */
	  double (*eta)(double, double); /** hard or soft thresholding */
	  double* (*sigextfct)(double*, int, int); /** extension fct: 'sym', 'zeros' */
} DenoisingParameters;

typedef struct{
	double theta1;
	double theta2;
} TimewarpParameters;

typedef struct{
	int *upath; /** contains y-coordinate in the path through the JxK matrix, length(upath)=K+J */
	int *spath;/** contains x-coordinate in the path through the JxK matrix */
	int J; /** length(u) */
	int K; /** length(s) */
} WarpPath;
	

/** \ingroup denoising */
/** Filling the fields is optional */
typedef struct {
	DenoisingParameters *den_params;
	TimewarpParameters  *tw_params;
	int N; /** num trials */
	int n; /** num points in each trial */
	double R; /** mean reaction time */
	double *Ri; /** individual reaction times in real time */
	double **sRi; /** individual reaction times in sampling points */
	double *times; /** sampling points -> real time */
	double **si; /** data */
	double *u; /** real signal*/
	double **ui; /** cleaned signals */
	void *additional; /** user data */
} ModelData;

#endif
