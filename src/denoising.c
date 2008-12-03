/* 
** denoising.c
**
** Made by (Matthias Ihrke)
** Login   <mihrke@localhost>
** 
** Started on  Thu Jul  5 23:30:56 2007 Matthias Ihrke
** Last update Sun Oct 14 12:31:45 2007 Matthias Ihrke
**/

#include "denoising.h"

/* ---------------------------------------------------------------------------- 
   -- Denoising routines                                                     -- 
   ---------------------------------------------------------------------------- */
/** Simple Running Median filter.
 * Is applied to each channel separately.
 */
void eeg_filter_running_median(EEGdata *s, int win){
  int i;
  for(i=0; i<s->nbchan; i++){
    running_median(s->d[i], s->n, win);
  }
}

/** Simple Running Median filter.
 * \code d[i] = median{ d[i-win], ..., d[i+win] } \endcode
 * -- choice of time-window 'win' is crucial
 * use 2*l+1 where l is the size of to-be-preserved details in sampling units
 * Ref: Schelter et al, 2005, chap. 6 (6.2.1)
 */
double* running_median(double *d, int n, int win){
  double *tmp, *dptr, *org;
  int i, awin, ewin;
  double med;

  org = (double*)malloc(n*sizeof(double));
  org = (double*)memcpy(org, d, n*sizeof(double));
  tmp = (double*)malloc((2*win+1)*sizeof(double));

  for(i=0; i<n; i++){
	 awin = i-win;
	 ewin = i+win;
	 if(awin<0) awin=0;
	 if(ewin>n-1) ewin=n;
    
    tmp = (double*)memcpy(tmp, &(org[awin]), sizeof(double)*(ewin-awin));
   
    gsl_sort(tmp, 1, ewin-awin);
    med = gsl_stats_median_from_sorted_data(tmp, 1, ewin-awin);

    d[i] = med;
  }
  free(tmp);
  free(org);
  return d;
}

/** Simple Moving Average filter in the time-domain.
 * s[i] = mean{ s[i-win], ..., s[i+win] }
 */
double* moving_average(double *s, int n, int win){
  int i, awin, ewin;
  double m;
  double *sptr, *org;

  org = (double*)malloc(n*sizeof(double));
  org = (double*)memcpy(org, s, n*sizeof(double));

  for(i=0; i<n; i++){
	 awin = i-win;
	 ewin = i+win;
	 if(awin<0) awin=0;
	 if(ewin>n-1) ewin=n;
	 m = gsl_stats_mean(&(org[awin]), 1, ewin-awin);
	 s[i] = m;
  }
  free(org);
  return s;
}

/**  Weighted Running Median filter.
 * Is applied to each channel separately.
 */
void eeg_filter_weighted_running_median(EEGdata *s, int win){
  int i;
  for(i=0; i<s->nbchan; i++){
    weighted_running_median(s->d[i], s->n, win, dist_euclidean);
  }
}
/** Weighted Running Median filter.
 * Ref: Schelter et al, 2005, chap. 6 (6.2.1)
 */
double* weighted_running_median(double *d, int n, int win, 
				double(*dist)(double,double)){
  double *pptr, *wptr;
  int i, j;
  int awin, ewin;
  double med;
  
  pptr = (double*)malloc((2*win+1)*sizeof(double));
  wptr = (double*)malloc((2*win+1)*sizeof(double));

  for(i=0; i<n; i++){
	 awin = i-win;
	 ewin = i+win;
	 if(awin<0) awin=0;
	 if(ewin>n-1) ewin=n;

    memcpy(pptr, &(d[awin]), sizeof(double)*(ewin-awin));
    for(j=0; j<(ewin-awin); j++){
      wptr[j] = (*dist)(i, awin+j);
    }
/* 	 dprintf("i=%i, awin=%i, ewin=%i\n", i, awin, ewin); */
    med = weighted_median_from_unsorted(pptr, wptr, ewin-awin);
	 
    d[i] = med;
  }
  free(pptr);
  free(wptr);

  return d;
}

/** Weighted Median computation.
 * Formula:
 *  \f[ WM = y_{(k)}\qquad k = \max\left\{h:\sum_{i=h}^n w_{(i)} \ge \frac{1}{2}\sum_{i=1}^{n} w_i\right\} \f]
 * \param d - data
 * \param w - corresponding weights (>=0)
 */
double weighted_median_from_unsorted(const double *d, const double *w, int n){
  size_t *permut; 
  double ref=0.0, refh;
  int k;
  int i, h, idx;

  permut=(size_t*)malloc(n*sizeof(size_t));
  gsl_sort_index(permut, d, 1, n);  
/*   dprintf("i\tperm[i]\td[i]\td[permut[i]]\n"); */
/*   for(i=0; i<n; i++){ */
/* 	 dprintf("%i\t%i\t%.2f\t%.2f\n", i, permut[i], d[i], d[permut[i]]); */
/*   } */

  for(i=0; i<n; i++)
    ref+=w[i];
  ref /= 2.0;
  
  for(h=n-1; h>0; h--){
    refh = 0.0;
    for(i=h; i<n; i++){ /* works because w_i>=0 */
      refh += w[permut[i]];
    }
    if(refh>=ref){
      k=h;
      break;
    }
  }
  idx = permut[k];
  free(permut);
  return d[idx];
}

double dist_euclidean(double x, double y){
  return fabs(x-y);
}

/** Wavelet estimation of single trial ERP's using Wang et al.'s
 * (2007) technique. 
 *
 * Formula \f[ \lambda = \sigma \sqrt{2\log_e{n\log_2{n}}} \f].
*/
double translation_invariant_thresholding(const double *data, int n){
  dprintf("Db: translation_invariant_thresholding\n");
  double sigma, lambda; /* population sd, threshold */
  
  sigma  = mad(data, n)/0.6745; /* funny constant, eh? */
  lambda = sigma * sqrt(2*log(n*glog(n, 2)));
  return lambda;
}

/** Wavelet estimation of single trial ERP's using the technique
 * labelled "conventional" in Wang et al. (2007).
 *
 * Formula \f[ \lambda = \sigma \sqrt{2\log_e{n}} \f].
 */
double conventional_thresholding(const double *data, int n){
  dprintf("Db: conventional_thresholding\n");
  double sigma, lambda; /* population sd, threshold */
  
  sigma  = mad(data,n)/0.6745; /* funny constant, eh? */
  lambda = sigma * sqrt(2*log(n));
  return lambda;
}
/** thresholding as proposed in the Matlab Wavelet-Toolbox.
    \code
       case 'heursure' 
         hthr = sqrt(2*log(n)); 
         eta = (norm(x).^2-n)/n; 
         crit = (log(n)/log(2))^(1.5)/sqrt(n); 
         if eta < crit 
             thr = hthr; 
         else 
             thr = min(thselect(x,'rigrsure'),hthr); 
         end 
	 \endcode
*/
double heuristic_sure(const double *data, int n){
  dprintf("Db: heuristic_sure\n");
  double hthr, etat, crit;
  hthr = sqrt(2*log(n));
  etat = (pow(vnorm(data, n, 2), 2)-n)/(double)n;
  crit = pow((log(n)/log(2)), 1.5)/sqrt(n);
  if(etat < crit)
    return hthr;
  else 
    return fmin(sureshrink(data, n), hthr);
}

/** Wavelet estimation of single trial ERP's using the procedure described
 * in Donoho et al., 1995. The implementation described in the outline
 * is used; */
double sureshrink(const double *data, int n){
  int i,k;
  double lambda, sigma, sure, suremin;
  double *tmp;
  dprintf("Db: sureshrink\n");

  tmp = (double*)malloc(n*sizeof(double));
  sigma  = mad(data, n)/0.6745;   

  for(i=0; i<n; i++)
    tmp[i] = fabs(data[i])/sigma;
  /*  tmp = memcpy(tmp, data, n*sizeof(double)); */
 
  /* compute the SURESHRINK threshold */
  suremin = DBL_MAX;
  qsort(tmp, n, sizeof(double), abscmp);
  
  for(k=0; k<n; k++){
    sure = n - 2*(k+1)+(n-k)*pow(fabs(tmp[k]), 2);
    for(i=0; i<k; i++)
      sure = sure + pow(fabs(tmp[i]), 2);
    if(sure<suremin){
      suremin = sure;
      lambda = fabs(tmp[k]);
    }
  }
  lambda = sigma * lambda;
  
  free(tmp);
  return lambda;
}

/** Generic Denoising: compute DWT of signal, call the thresholding
   function 'threshfct' for each resolution level and IDWT the signal */
int generic_denoising(double *data, int n, int L, 
		      double(*threshfct)(const double*, int), double(*etafct)(double,double)){
  gsl_wavelet *w;
  gsl_wavelet_workspace *work;
  int j, J, k, offset;
  double lambda; /* population sd, threshold */
  dprintf("Db: generic_denoising\n");

  w = gsl_wavelet_alloc(gsl_wavelet_daubechies, 20);
  work = gsl_wavelet_workspace_alloc(n);

  gsl_wavelet_transform_forward(w, data, 1, n, work);

  /* -- thresholding here -- */
  J = (int)round(glog((double)n, 2));
  for(j=L; j<J; j++){ /* loop through levels */
    offset = (int)pow(2, j);
    lambda = (*threshfct)(&(data[offset]), offset);

    for(k=offset; k<2*offset; k++){ /* loop through coefficients */
      /* soft or hard thresholding */
      data[k]=(*etafct)(data[k], lambda);
    }
  }
  /* -- thresholding end -- */

  gsl_wavelet_transform_inverse(w, data, 1, n, work);
  gsl_wavelet_free(w);
  gsl_wavelet_workspace_free(work);
  return 0;
}

/** Extend data to length 2^j using sigextfct and denoise it.
 * \see generic_denoising()
 */
int extend_and_denoise(double *data, int n, int L, 
		       double(*threshfct)(const double*, int), 
		       double(*etafct)(double,double), 
		       double*(*sigextfct)(double*, int, int)){
  int J;
  double *tmp, *dptr;

  J = (int)round(glog((double)n, 2));
  J++;
  if(((int)pow(2, J))<n)
    dprintf(" Warning: Signal length %i is not a power of 2 -> extending signal to length %i\n",
	    n, (int)pow(2,J));
  tmp = (double*)malloc((int)pow(2,J) * sizeof(double));
  tmp = memcpy(tmp, data, n*sizeof(double));
  dptr=(*sigextfct)(tmp, n, (int)pow(2, J));
  generic_denoising(tmp, (int)pow(2,J), L, threshfct, etafct);
  
  data = memcpy(data, dptr, n*sizeof(double));

  free(tmp);
  return 0;
}

/** Implements soft-thresholding.
 * Formula: \f$ \eta_s(\lambda, w)=\f$
 * \ingroup thresholding
 */
double eta_s(double d, double lambda){
  /* soft thresholding, formula (6) */
/*   dprintf("Db: eta_s\n"); */
  if(fabs(d)<=lambda){
    d=0.0;
  } else if(d>lambda){
    d=d-lambda;
  } else if(d<(-1)*lambda){
    d=d+lambda;
  }
  return d;
}
/** Implements hard-thresholding.
 * Formula: \f$ \eta_h(\lambda, w)=\f$
 * \ingroup thresholding
 */
double eta_h(double d, double lambda){
  /* hard thresholding */
/*   dprintf("Db: eta_h\n"); */
  if(fabs(d)<=lambda)
    d=0.0;
  return d;
}


/** Extend data to length 2^j using sigextfct and denoise it.
 * \see generic_denoising()
 * \see extend_and_denoise()
 * data is directly written into the eegdata-struct
 */
void eeg_wavelet_denoise(EEGdata *eeg, int L, 
								double(*threshfct)(const double*, int), 
								double(*etafct)(double,double), 
								double*(*sigextfct)(double*, int, int)){
  int c;
  
  for( c=0; c<eeg->nbchan; c++ ){
	 extend_and_denoise( eeg->d[c], eeg->n, L, threshfct, etafct, sigextfct );
  }
}
