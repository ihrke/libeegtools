/***************************************************************************
 *   Copyright (C) 2008 by Matthias Ihrke   *
 *   mihrke@uni-goettingen.de   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "mathadd.h"

double glog(double v, int b){
  /* compute log_b(v) using ansi-C log */
  return (double)((double)log((double)v)/(double)log((double)b));
}

double mad(const double *data, int n){
  /* For a univariate data set X = X_1, X_2, ..., X_n, the MAD is defined as 
	MAD = median(|X_i - median(X)|)
	(this function has been validated against Matlab's routine, i.e. it works)
  */
	double *tmp;
	double med;
	int i;
	tmp = (double*)malloc(n*sizeof(double));
	tmp = (double*)memcpy((void*)tmp, (const void*)data, sizeof(double)*n);

	gsl_sort(tmp, 1, n);
	med = gsl_stats_median_from_sorted_data(tmp, 1, n);
	for(i=0; i<n; i++)
		tmp[i] = fabs(data[i]-med);
	gsl_sort(tmp, 1, n);
	med = gsl_stats_median_from_sorted_data(tmp, 1, n);
  
	free(tmp);
	return med;
}


int    abscmp(const void *p1, const void *p2){
	if(fabs(*(double*)p1)<fabs(*(double*)p2))
		return -1;
	else if(fabs(*(double*)p1)>fabs(*(double*)p2))
		return 1;
	else 
		return 0;
}

double vnorm(const double *v, int n, int p){
	/* vector norm (p) */
	int i;
	double norm = 0;
	for(i=0; i<n; i++)
		norm = norm + pow(fabs(v[i]), p);
	return pow(norm, 1/(double)p);
}


double maxel(double *v, int n){
  int i;
  double m=DBL_MIN;
  for(i=0; i<n; i++)
    if(v[i]>m) m=v[i];
  return m;
}

int maxeli(int *v, int n){
  int i;
  int m=INT_MIN;
  for(i=0; i<n; i++)
    if(v[i]<m) m=v[i];
  return m;
}

int closest_index(double *v, int n, double c){
  int i, index=0;
  double tmp=DBL_MAX;
  for(i=0; i<n; i++)
    if(fabs(v[i]-c) < tmp){
      index=i;
      tmp = fabs(v[i]-c);
    }
  return index;
}

double* sampled_line(double *ntimes, int n, double start, double end){
  /* return pointer to *ntimes with n equally spaced doubles, where
     ntimes[0]=start and ntimes[n-1]=end */
  double step;
  int i;

  step=(end-start)/(double)(n-1);
  for(i=0; i<n; i++)
    ntimes[i]=start+(i*step);

  return ntimes;
}

/** same as first:last in matlab.
 * Data ist allocated in the function and needs to be freed
 * by the caller.
 */
int* linspace(int first, int last){
	int i;
	int *s;
	massert(first>last, "first>last: %i>%i", first, last);
	s=(int*)malloc((last-first+1)*sizeof(int));
	int j=0;
	for(i=first; i<last; i++){
		s[j]=i;
		j++;
	}
	return s;
}

/** performs linear interpolation for the values x2, given
 * x1 and y1. A pointer to the interpolated values y2 is 
 * returned. The x2 must ly within the x1, the x1 must be
 * increasing in magnitude. 
 */
double* lininterp(const double *x1, const double *y1, int n1, 
		  const double *x2,       double *y2, int n2){
  int i, j;
  double m, n;

  massert(x2[0]<x1[0], "x2[0]<x1[0]: %f < %f", x2[0], x1[0]);
  massert(x2[n2-1]>x1[n1-1], "x2[n2-1]>x1[n1-1]: %f > %f", x2[n2-1], x1[n1-1]);
/*   fprintf(stderr, "n1=%i, n2=%i\n", n1, n2); */
/*   vprint_vector("x1", x1, n1); */
/*   vprint_vector("y1", y1, n1); */
/*   vprint_vector("x2", x2, n2); */
  for(i=0; i<n2; i++){
    j=0; while(j<n1 && x2[i]>x1[j]) j++; /* set j */
/*     printf("i=%i (%2.2f), j=%i (%2.2f)\n",i, x2[i], j, x1[j]); */
    if(j==n1) n--;
    m = (y1[j]-y1[j-1])/(x1[j]-x1[j-1]);
    n = y1[j]-m*x1[j];
    y2[i] = m*x2[i]+n;
  }
  return y2;
}
/** Leave one out cross validation.
 * \param u reference signal
 * \param si the data
 * \param err the estimate for the external prediction error (curve); must be allocated by caller
 * \param N num trials (si)
 * \param n num data points in u and each si
 * \param model the function to estimate u from si. A model takes the si, an empty array with n positions to put the estimate, an integer array giving the indices in the si (to match them to the userdata, kappa), n and userdata.
 * \param userdata is simply passed to the model (could be a struct with everything needed)
 */
double* loocv(const ModelData *m, double* err, double*(*model)(const ModelData*,double*)){
	double **newsi, *newRi;
	ModelData *mtmp;
	double *hatu;
	int i,j, k, n, N;
	n = m->n;
	N = m->N;
	newsi = (double**)malloc(2*N*sizeof(double*));
	hatu = (double*)malloc(n*sizeof(double));
	newRi = (double*)malloc(2*N*sizeof(double));
	mtmp = (ModelData*)malloc(sizeof(ModelData));
	
	dprintf("n1=%i\n", m->n);
	for(i=0; i<n; i++) err[i]=0.0;
	for(i=0; i<2*N; i++){
		newsi[i]=m->si[i%N];
		newRi[i]=m->Ri[i%N];
	}
	for(i=0; i<N; i++){
		copy_modeldata( m, mtmp );
		mtmp->N = N-1;
		mtmp->Ri = &(newRi);
		mtmp->si = &(newsi[i]);
		hatu = (*model)(mtmp, hatu);
		for(j=0; j<n; j++)
			err[j] += SQR((m->u)[j]-hatu[j]);
	}
	for(i=0; i<n; i++) err[i]/=N;
	
	free(newsi);
	free(newRi);
	free(hatu);
	return err;
}	 
