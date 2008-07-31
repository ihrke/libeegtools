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
 * - writer.c writer.h  - write data
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
 * - Some demo programs are provided in ./progs of the distribution. 
 * - CAUTION: some of the functions are not well tested (alpha). 
 *
 * \section Compilation
 * - You need the GNU Scientific Library >=1.09, available from http://www.gnu.org/software/gsl/
 * - then compilation and installation according to GNU autotools:
\code
$ tar xvfz libeegtools-ver.tar.gz
$ cd libeegtools-ver.tar.gz
$ ./configure
$ make
$ make install
\endcode
 * 
 * - to install the matlab-files as well, do
\code
$ cd matlab
$ make -f Makefile.matlab
\endcode
 * this will compile the Matlab-Mex-files. You can then use all MATLAB-functions in 
 * this directory from the MATLAB-command line.
 * - after this, you can install the eeglab-plugin by copying all files from ./matlab and 
 *   ./eeglab_plugin to a directory in your eeglab-plugin path.
 *
 * 
\bf History:
\code
 * $Log$
 * Revision 1.7  2008/07/31 20:07:53  mihrke
 * evidence for different warping :-)
 *
 * Revision 1.5  2008/05/30 19:53:09  mihrke
 * version 0.2a
 *
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

#ifdef __cplusplus
extern "C" {
#endif

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
	 int *upath; /** contains y-coordinate in the path through the JxK matrix, maximal length(upath)=K+J */
	 int *spath;/** contains x-coordinate in the path through the JxK matrix */
	 int J; /** length(u) (not length(upath)) */
	 int K; /** length(s) (not length(spath)) */
	 //	 int real_length; /** real length of the path (maximal K+J) */
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


  typedef struct{
	 char label[255];
	 double x;
	 double y;
	 double z;
  } ChannelInfo;

  static ChannelInfo channelcoords_64[] = {
	 {"Fp1",80.8,26.1,-4},
	 {"Fp2",80.8,-26.1,-4},
	 {"F3",57.6,48.2,39.9},
	 {"F4",57.6,-48.1,39.9},
	 {"C3",3.87e-15,63.2,56.9},
	 {"C4",3.87e-15,-63.2,56.9},
	 {"P3",-57.6,48.2,39.9},
	 {"P4",-57.6,-48.1,39.9},
	 {"EOGru",100,-40,-80},
	 {"EOGro",100,-70,-40},
	 {"F7",49.9,68.4,-7.49},
	 {"F8",49.9,-68.4,-7.49},
	 {"T7",5.18e-15,84.5,-8.85},
	 {"T8",5.18e-15,-84.5,-8.85},
	 {"P7",-49.9,68.4,-7.49},
	 {"P8",-49.9,-68.4,-7.49},
	 {"Fz",60.7,0,59.5},
	 {"Cz",5.2e-15,0,85},
	 {"Pz",-60.7,-7.44e-15,59.5},
	 {"Oz",-85,-1.04e-14,-1.79},
	 {"FC1",32.4,32.4,71.6},
	 {"FC2",32.4,-32.4,71.6},
	 {"CP1",-32.4,32.4,71.6},
	 {"CP2",-32.4,-32.4,71.6},
	 {"FC5",28.8,76.2,24.2},
	 {"FC6",28.8,-76.2,24.2},
	 {"CP5",-28.8,76.2,24.2},
	 {"CP6",-28.8,-76.2,24.2},
	 {"TP9",-23.3,70.1,-42.1},
	 {"TP10",-23.3,-70.1,-42.1},
	 {"EOGlu",100,40,-80},
	 {"ECGlo",100,70,-40},
	 {"F1",59.9,26,54.4},
	 {"F2",59.9,-26,54.4},
	 {"C1",2.11e-15,34.5,77.7},
	 {"C2",2.12e-15,-34.6,77.6},
	 {"P1",-59.9,26,54.4},
	 {"P2",-59.9,-26,54.4},
	 {"AF3",76.2,31.5,20.8},
	 {"AF4",76.2,-31.5,20.8},
	 {"FC3",31,59.3,52.5},
	 {"FC4",31,-59.3,52.5},
	 {"CP3",-31,59.3,52.5},
	 {"CP4",-31,-59.3,52.5},
	 {"PO3",-76.2,31.5,20.8},
	 {"PO4",-76.2,-31.5,20.8},
	 {"F5",54,63.1,18.1},
	 {"F6",54,-63,18.2},
	 {"C5",4.95e-15,80.8,26.3},
	 {"C6",4.95e-15,-80.8,26.3},
	 {"P5",-54,63.1,18.1},
	 {"P6",-54,-63,18.2},
	 {"AF7",68.7,49.7,-5.96},
	 {"AF8",68.7,-49.7,-5.95},
	 {"FT7",26.2,80.4,-8.51},
	 {"FT8",26.2,-80.4,-8.51},
	 {"TP7",-26.2,80.4,-8.51},
	 {"TP8",-26.3,-80.4,-8.51},
	 {"PO7",-68.7,49.7,-5.96},
	 {"PO8",-68.7,-49.7,-5.95},
	 {"Fpz",85,0,-1.79},
	 {"AFz",79,0,31.3},
	 {"CPz",-32.9,-4.03e-15,78.4},
	 {"POz",-79,-9.68e-15,31.3}
  };


#ifdef __cplusplus
}
#endif

#endif
