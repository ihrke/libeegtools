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
 * d[i] = median{ d[i-win], ..., d[i+win] }
 * -- choice of time-window 'win' is crucial
 * use 2*l+1 where l is the size of to-be-preserved details in sampling units
 * Ref: Schelter et al, 2005, chap. 6 (6.2.1)
 */
double* running_median(double *d, int n, int win){
  double *tmp, *dptr;
  int i, cwin;
  double med;

  if(win % 2==0) win++;
  tmp = (double*)malloc(win*sizeof(double));

  for(i=0; i<n; i++){
    dptr = &(d[i]);
    cwin = win;
    if(i-win/2 < 0){
      dptr=d;
      cwin = i+win/2;
    } 
    if(i+win/2 >= n) cwin=n-i;
    
    tmp = (double*)memcpy((void*)tmp, (const void*)(dptr), sizeof(double)*cwin);
   
    gsl_sort(tmp, 1, cwin);
    med = gsl_stats_median_from_sorted_data(tmp, 1, cwin);

    d[i] = med;
  }
  free(tmp);
  return d;
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
  double *pptr, *wptr, *dptr;
  int i, j;
  int cwin, cstart;
  double med;
  
  if(win % 2==0) win++;
  pptr = (double*)malloc(win*sizeof(double));
  wptr = (double*)malloc(win*sizeof(double));

  for(i=0; i<n; i++){
    dptr = &(d[i]);
    cwin = win;
    cstart = i-win/2;
    if(i-win/2 < 0){
      dptr=d;
      cstart = 0;
      cwin = i+win/2;
    } 
    if(i+win/2 >= n) cwin=n-i;
    
    pptr = (double*)memcpy((void*)pptr, (const void*)(dptr), sizeof(double)*cwin);
    for(j=0; j<cwin; j++){
      wptr[j] = (*dist)(i, cstart+j);
    }
   
    med = weighted_median_from_unsorted(pptr, wptr, cwin);

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
  int *permut; 
  double ref=0.0, refh;
  int k;
  int i, h;

  permut=(int*)malloc(n*sizeof(int));
  gsl_sort_index(permut, d, 1, n);  
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
  free(permut);

  return d[permut[k]];
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

/*   char tmpstring[255]; */
/*   extern Engine *matlab;  */
/*   ml_plot(matlab, NULL, data, n, "r", 1);  */
  
  /* -- thresholding here -- */
  J = (int)round(glog((double)n, 2));
  for(j=L; j<J; j++){ /* loop through levels */
    offset = (int)pow(2, j);
    lambda = (*threshfct)(&(data[offset]), offset);

/*     fprintf(stderr, " j=%i, range %f - %i\n", j, pow(2,(double)(j)),(int)pow(2,j+1)); */
/*     sprintf(tmpstring, "hold on; plot(%i:%i, %f, 'g');",(int)pow(2,j),(int)pow(2,j+1), lambda); */
/*     engEvalString(matlab, tmpstring); */

    for(k=offset; k<2*offset; k++){ /* loop through coefficients */
      /* soft or hard thresholding */
      data[k]=(*etafct)(data[k], lambda);
    }
  }
/*   ml_wait(matlab); */
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

double eta_h(double d, double lambda){
  /* hard thresholding */
/*   dprintf("Db: eta_h\n"); */
  if(fabs(d)<=lambda)
    d=0.0;
  return d;
}

/* ---------------------------------------------------------------------------- 
   -- Merit Measures                                                         -- 
   ---------------------------------------------------------------------------- */

/** \f$ 
    RMSE = \sqrt{ \frac{1}{\#r}  \sum{ (r-d)^2)}} 
    \f$ 
 */
double rmse(const double *r, const double *d, int n){
  int i;
  double res=0.0;
  for(i=0; i<n; i++)
    res = res + pow(r[i]-d[i], 2);
  return sqrt(res*1/n);
}

/** \f$
    SNR = 10\log_{10} \frac{\sum{r^2}}{\sum{(r-d)^2}} 
 \f$*/
double snr (const double *r, const double *d, int n){
  int i;
  double res, tmp;
  res = 0.0; tmp = 0.0;
  for(i=0; i<n; i++){
    res = res+pow(r[i], 2);
    tmp = tmp+pow(r[i]-d[i], 2);
  }
  return 10*(glog(res/tmp, 10));
}

/* ---------------------------------------------------------------------------- 
   -- Signal extension routines                                              -- 
   ---------------------------------------------------------------------------- */
/* the extension functions return a pointer to the former data[0],
   because this is where the unextended signal began;
   Example:
      sigext([1 2 3 - - - -]) -> [0 0 1 2 3 0 0] 
                                      ^ ptr
   Assumptions (not for full generality!):
      1) ns <= n
      2) n <= 2*ns
*/

/** \code [1 2 3 - - -] -> [0 1 2 3 0 0] \endcode*/
double* sigext_zeros(double *data, int ns, int n){
  int offset, i=0; /* for signal */
  double *dptr;
  dprintf("Db: sigext_zeros\n");

  offset = (n-ns)/2;
  dptr = memmove(&(data[offset]), data, ns*sizeof(double));
  for(i=1; i<=offset; i++){ /*( ((n-ns)%2==0) ? offset : offset-1); i++){*/
    data[offset-i]=0.0; 
    data[offset+ns+i-1]=0.0;
  }
  data[n-1]=0.0;
  return dptr;
}


/** \code [1 2 3 - - -] -> [1 2 3 0 0 0] \endcode */
double* sigext_zerosr(double *data, int ns, int n){
  dprintf("Db: sigext_zerosr\n");
  int i; /* for signal */
  for(i=ns; i<n; i++) data[i]=0.0; 
  return data;
}

/** \code [1 2 3 - - - - -] -> [2 1 1 2 3 3 2 1] \endcode*/
double* sigext_sym(double *data, int ns, int n){
  int offset, i=0; 
  double *dptr;
  dprintf("Db: sigext_sym\n");

  offset = (n-ns)/2;
  dptr = memmove(&(data[offset]), data, ns*sizeof(double));
  for(i=1; i<=offset; i++){ /*( ((n-ns)%2==0) ? offset : offset-1); i++){*/
    data[offset-i]=data[offset+i-1]; 
    data[offset+ns+i-1]=data[offset+ns-i];
  }
  data[n-1]=((ns-n)%2==0 ? data[offset+(n-ns)/2] : data[offset+(n-ns)/2-1]);
  return dptr;
  
}

double* sigext_smooth(double *data, int ns, int n){
  /* [1 2 3 - - - - -] -> [1 1 1 2 3 3 3 3] */
  int offset, i=0; /* for signal */
  double *dptr;
  dprintf("Db: sigext_smooth\n");

  offset = (n-ns)/2;
  dptr = memmove(&(data[offset]), data, ns*sizeof(double));
  for(i=1; i<=offset; i++){ /*( ((n-ns)%2==0) ? offset : offset-1); i++){*/
    data[offset-i]=data[offset]; 
    data[offset+ns+i-1]=data[ns-1];
  }
  data[n-1]=data[n-2];
  return dptr;
}
