/** \mainpage 
 *
 * This is libEEGTools.
 *
 * Synopsis:
 * - helper.h, helper.c contain helper math, print, plot functions etc.
 * - denoising.h, denoising.c - core functionality
 *
 * History:
\code
 * $Log$
 * Revision 1.3  2008/03/14 15:44:27  mihrke
 * distribution 0.15
 *
\endcode
 */

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

typedef struct{
  int nbchan; /** number of channels */
  int n; /** number of samples */
  double **d; /** data */
} EEGdata;


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
