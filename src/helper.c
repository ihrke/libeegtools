/* helper.c
 *
 */	
#include "helper.h"

/* ---------------------------------------------------------------------------- 
   -- Helper functions                                                       -- 
   ---------------------------------------------------------------------------- */

/** deep copy of double ptrptr
 */
double** copy_double_ptrptr(const double **s, int N, int n){
  double **r;
  int i, j;

  r = (double**)malloc(N*sizeof(double*));
  for(i=0; i<N; i++){
    r[i] = (double*)malloc(n*sizeof(double));
	 for(j=0; j<n; j++){
		r[i][j] = s[i][j];
	 }
  }
  return r;
}

/** shallow copy of ModelData struct
 */
void    copy_modeldata(const ModelData *m1, ModelData *m2){
	m2->den_params=m1->den_params;
	m2->tw_params=m1->tw_params;
	m2->N=m1->N;
	m2->n=m1->n;
	m2->R=m1->R;
	m2->Ri=m1->Ri;
	m2->sRi=m1->sRi;
	m2->times=m1->times;
	m2->si=m1->si;
	m2->u=m1->u;
	m2->ui=m1->ui;
	m2->additional=m1->additional;
}

void free_warppath(WarpPath *p){
	free(p->upath);
	free(p->spath);
	free(p);
}

void print_modeldata(FILE *out, const ModelData *m){
	fprintf(out, "ModelData struct:\n");
	fprintf(out, "  den_params=\n");
	print_denoisingparameters(out, m->den_params);
	fprintf(out, "  tw_params=\n");
	print_timewarpparameters(out, m->tw_params);
	fprintf(out, "  N=%i\n", m->N);
	fprintf(out, "  n=%i\n", m->n);
	fprintf(out, "  R=%i\n", m->R);
	fprintf(out, "  Ri=%p\n", m->Ri);
	fprintf(out, "  sRi=%p\n", m->sRi);
	fprintf(out, "  times=%p\n", m->times);
	fprintf(out, "  si=%p\n", m->si);
	fprintf(out, "  u=%p\n", m->u);
	fprintf(out, "  ui=%p\n", m->ui);
	fprintf(out, "  additional=%p\n", m->additional);
}
void    print_denoisingparameters(FILE *out, const DenoisingParameters *p){
	fprintf(out, "   DenoisingParameters struct:\n"	);
	fprintf(out, "      L=%i\n", p->L);
	fprintf(out, "      sigextfct=%p\n", p->sigextfct);
	fprintf(out, "      cleanfct=%p\n", p->cleanfct);
	fprintf(out, "      eta=%p\n", p->eta);
}
void    print_timewarpparameters(FILE *out, const TimewarpParameters *p){
	fprintf(out, "   TimewarpParameters struct:\n"	);
	fprintf(out, "      theta1=%i\n", p->theta1);
	fprintf(out, "      theta1=%i\n", p->theta2);
}

		
void free_modeldata(ModelData *m){
	int i;
	if(!(m->den_params==NULL))	free(m->den_params);
	if(!(m->tw_params==NULL))	free(m->tw_params);
	if(!(m->additional==NULL))	free(m->additional);
	if(!(m->times==NULL)) free(m->times);
	if(!(m->sRi==NULL)) free(m->sRi);
	if(!(m->Ri==NULL)) free(m->Ri);
	if(!(m->u==NULL)) free(m->u);
	for(i=0; i<m->N; i++){
		if(m->si!=NULL)
			if(!(m->si[i]==NULL)) free(m->si[i]);
		if(m->ui!=NULL)
			if(!(m->ui[i]==NULL)) free(m->ui[i]);
	}
	if(m->si!=NULL) free(m->si);
	if(m->ui!=NULL) free(m->ui);
	free(m);
}

int v_printf(int v, char *format, ...){
  extern int *verbosity;
  va_list arglist;

  va_start(arglist,format);

  if (v<=*verbosity)
    vfprintf(stderr,format,arglist);
  va_end(arglist);

  return 0;
}

int vprint_vector(const char* name, double *v, int n){
  int i;
  fprintf(stderr, "%s = \n", name);
  for(i=0; i<n; i++)
    fprintf(stderr, " %2.2f ", v[i]);
  fprintf(stderr, "\n");
  return 0;
}


void   errormsg(int err_no, int fatal){
  switch(err_no){
  case ERR_IO:
    v_printf(0, "IO Error");
    break;
  case ERR_GSL:
    v_printf(0, "Error in the GSL-library");
    break;
  case ERR_PLOT:
    v_printf(0, "Error in the Plot-library");
    break;
  default:
    v_printf(0, "Unknown Error number");
  }
  if(fatal){
    v_printf(0, "... Fatal\n");
    exit(err_no);
  } else v_printf(0, "\n");
}

#ifdef HAVE_PLOTLIB
void plot_coordsys(double xmin, double xmax, double ymin, double ymax, double *labels, int numlab){
  int i;
  char tmp[255];
  pl_pencolorname("black");
  pl_linewidth(0.1);
  pl_fline(xmin, 0.0, xmax, 0.0);
  pl_fline(0.0, YSCALE*ymin, 0.0, YSCALE*ymax);
  pl_fontname("Courier");
  pl_ffontsize(30);
  for(i=0; i<numlab; i++){
    pl_fmove(labels[i], 0.0);
    sprintf(tmp, "%.0f", labels[i]);
    pl_alabel('c', 't', tmp);
  }
}

int plot_init(double xmin, double xmax, double ymin, double ymax, char* plotdev){
  int h, i, numlabels;
  double *labels;
  char tmp[20];
  if(!plotdev)
    return -1;
  sprintf(tmp, "%ix%i", SCREENX, SCREENY);
  
  numlabels = (xmax-xmin)/100.0;
  labels = (double*)malloc(numlabels*sizeof(double));
  labels[0]=xmin;
  for(i=1; i<numlabels; i++)
    labels[i]=labels[i-1]+100.0;

  pl_parampl("BITMAPSIZE", tmp);
  pl_parampl("PAGESIZE", "a4");
  if((h=pl_newpl(plotdev, stdin, stdout, stderr))<0){
    errormsg(ERR_PLOT, 0);
    return ERR_PLOT;
  }
  pl_selectpl(h);
  if(pl_openpl()<0){
    errormsg(ERR_PLOT,0);
    return ERR_PLOT;
  }
  
  //pl_space(xmin, ymin*SCREENX/SCREENY, xmax, ymax*SCREENX/SCREENY);
  pl_space(xmin, xmin, xmax, xmax);
  plot_coordsys(xmin, xmax, ymin, ymax, labels, numlabels); 

  return h;
}

void   plot_trace(const double *times, const double *d, int n, const char *color){
  /* if times==NULL then the points are plotted from zero on;
   * else both arrays must be of length n
   */
  int i;
  pl_endpath();
  pl_pencolorname(color);
  pl_linewidth(0.1);
  if(!times)
    pl_fmove(0.0, d[0]);
  else
    pl_fmove(times[0], d[0]);

  for(i=0; i<n; i++){
    if(!times)
      pl_fcont((double)i, d[i]); 
    else
      pl_fcont(times[i], YSCALE*d[i]);
  }
  pl_endpath();
}

void plot_close(int handle){
  if(pl_closepl()<0)
    errormsg(ERR_PLOT, 0);
  pl_selectpl(0);
  if(pl_deletepl(handle)<0)
    errormsg(ERR_PLOT, 0);
}
#endif



#ifdef HAVE_MATLAB
Engine* ml_init(){
  Engine *m;
  if (!(m = engOpen(MATLAB_STARTUP_CMD))) {
    dprintf("\nCan't start MATLAB engine\n");
    return NULL;
  }
  dprintf("Matlab successfully started...\n");
  return m;
}

int ml_close(Engine *m){
  return engClose(m);
}


void ml_plot(Engine *matlab, const double *r, const double *v, int n, const char *color, int new){
  mxArray *a1, *a2;
  double *d1,*d2;
  char tmpstring[255];
 
  if(r!=NULL){
    a1 = mxCreateDoubleMatrix(1,n,mxREAL);
    d1 = mxGetPr(a1); 
    d1 = memcpy(d1, r, n*sizeof(double)); 
    if(engPutVariable(matlab, "times", a1))
       fprintf(stderr,"Error putting variable to MATLAB\n");
  }
  a2 = mxCreateDoubleMatrix(1,n,mxREAL);
  d2 = mxGetPr(a2); d2 = memcpy(d2, v, n*sizeof(double));

  if(engPutVariable(matlab, "sig", a2))
    fprintf(stderr,"Error putting variable to MATLAB\n");
  
  if(new) engEvalString(matlab, "figure;");
  if(r!=NULL)
    sprintf(tmpstring, "plot(times, sig, '%s');", color);
  else
    sprintf(tmpstring, "plot(sig, '%s');", color);

  dprintf("Evaluating: '%s'\n",tmpstring); 
  engEvalString(matlab, tmpstring);

  if(r!=NULL)  mxDestroyArray(a1);
  mxDestroyArray(a2);
  return;
} 

void ml_plot_path(Engine *matlab, const int *path, int K){
  int J, m, i, j, k;
  double *x, *y;
  mxArray *xp, *yp;
  char tmpstring[255];

  J = path[K-1];
  m = MAX(J, K);

  xp = mxCreateDoubleMatrix(1,m,mxREAL);
  x = mxGetPr(xp); 

  yp = mxCreateDoubleMatrix(1,m,mxREAL);
  y = mxGetPr(yp); 

  i=0; j=0; k=0; 
  for(i=0; i<m; i++){
    x[i] = (double)k;
    y[i] = (double)j;

    if(path[k]==j){  k++;  }
    else if(path[k]>j) j++;
  }

  if(engPutVariable(matlab, "xpath", xp) ||engPutVariable(matlab, "ypath", yp))
    fprintf(stderr,"Error putting variable to MATLAB\n");
  sprintf(tmpstring, "figure; plot(xpath, ypath); axis([0 %i 0 %i]);", K, J);
  engEvalString(matlab,tmpstring);
  
  mxDestroyArray(xp);  mxDestroyArray(yp);
  return;
}

void ml_wait(Engine *matlab){
  engEvalString(matlab, "waitforbuttonpress;");
}
#endif
