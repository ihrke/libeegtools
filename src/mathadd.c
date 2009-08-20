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

/** return index in vector that is closest to a value
 */	  
int closest_index(const double *v, int n, double c){
  int i, index=-1;
  double tmp=DBL_MAX;
  /* dprintf( " v[%i] = %f, n=%i\n", i, v[i], n );	 */
  for(i=0; i<n; i++){
    if(fabs(v[i]-c) < tmp){
      index=i;
      tmp = fabs(v[i]-c);
    }
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
/** same as first:step:last in matlab.
 * A sequence first, first+step, first+(2*step), ... ,last
   is computed. The length of this  sequence is 
	\f$ n=\lfloor \frac{last-first}{step} \rfloor+1 \f$
	\param v either enough free allocated space, or ALLOC_IN_FCT (NULL) 
	         - data ist allocated in the function and needs to be freed  by the caller.
	\param n - number of doubles in v (return value)
 */
double* linspace_dbl(double first, double last, double step, double *v, int *n){
  int i;
	massert(first>last, "first>last: %i>%i", first, last);
	*n = (int) round(((last-first)/step)+1);
	if( v==ALLOC_IN_FCT ){
	  v=(double*)malloc((*n)*sizeof(double));
	}
	for( i=0; i<*n; i++ ){
	  v[i] = first+i*step;
	}

	return v;
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
	for(i=first; i<=last; i++){
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

/* ---------------------------------------------------------------------------- 
   -- vector ops                                                             -- 
   ---------------------------------------------------------------------------- */

/** find minimum element in vector.
	 \param v vector of length n
	 \param idx the index in v where v[i] is minimal
	 \return v[i]
 */
double vector_min( double *v, int n, int *idx ){
  int i;
  double min=DBL_MAX;
  *idx=-1;

  for( i=0; i<n; i++ ){
	 if( v[i] < min ){
		min = v[i];
		*idx = i;
	 }
  }

  return min;
}
/** find max element in vector.
	 \param v vector of length n
	 \param idx the index in v where v[i] is max
	 \return v[i]
 */
double vector_max( double *v, int n, int *idx ){
  int i;
  double max=DBL_MIN;
  *idx=-1;

  for( i=0; i<n; i++ ){
	 if( v[i] > max ){
		max = v[i];
		*idx = i;
	 }
  }

  return max;
}

/** create a random permutation of the elements in permut.
	 This is crude and simple, the function graps n pairs of indices and swaps them.
 */
void vector_shuffle_int( int *permut, int n ){
  int i1,i2, i;

  for( i=0; i<n; i++){
	 i1 = (random() / (RAND_MAX / n+1));
	 i2 = (random() / (RAND_MAX / n+1));
	 swap2i( &(permut[i1]), &(permut[i2]));
  }
}


/** remove a scalar from all values in a vector.
 */
void  vector_minus_scalar( double *v, int n, double val ){
  int i;
  for( i=0; i<n; i++ ){
	 v[i] = v[i]-val;
  }
}

/** \param v - if NULL, own memory is allocated
	 \param n - length(v)
	 \param val - set v[i]=val for all i=0,...,n-1
	 \return v
 */
double* vector_init( double *v, int n, double val ){
  int i;
  if( v==NULL){
	 v = (double*)malloc( n*sizeof(double) );
  }
  for(i=0; i<n; i++){
	 v[i] = val;
  }

  return v;
}

/** calculate \f$ |\vec{v_1} - \vec{v_2}|^2 \f$ 
 */
double vector_euclidean_distance( const double *v1, const double *v2, int n ){
  int i;
  double r=0;
  for( i=0; i<n; i++ ){
	 r += SQR( v1[i]-v2[i] );
  }
  return sqrt(r);
}

/* ---------------------------------------------------------------------------- 
   -- Matrix ops                                                             -- 
   ---------------------------------------------------------------------------- */

/** fills m with uniformely drawn random values from [lower, upper].
 */
double** matrix_rand( double **m, int N, int M, double lower, double upper ){
  int i,j;
  
  srand((long)time(NULL));
  if( !m ){
	 m = matrix_init( N, M );
  }
  for( i=0; i<N; i++ ){
	 for( j=0; j<M; j++ ){
		m[i][j] = ( (((double)rand()) / RAND_MAX)*(upper-lower))+lower;
	 }
  }

  return m;
}
/** Normalizes matrix by its maximum value:
	 \f[
	 \hat{M} = \frac{1}{\mbox{max}(M)} M
	 \f]
 */
void     matrix_normalize_by_max( double **m, int M, int N ){
  double max; 

  max = matrix_max( m, M, N, NULL, NULL );
  matrix_divide_scalar( m, M, N, max);
}

/** Delete a row in a matrix. Memory remains allocated and the row pointer
 *  is moved to the end of the matrix. Index runs from 0,...,N-1
 * \param m - matrix 
 * \param N - num of rows
 * \param n - num of cols
 * \param row- num of row to delete
 */
double** matrix_delrow(double **m, int N, int n, int row){
  double *tmp;
  int i;
  massert(row>N-1, "cannot delete row %i because only %i rows\n", row, N);
  tmp = m[row];
  for(i=row; i<N-1; i++){
    m[i]=m[i+1];
  }
  m[N-1] = tmp;
  return m;
}

/** return a 0-initialized matrix of dimension NxM
 */
double** matrix_init(int N, int M){
  int i,j;
  double **d;
  d = (double**) malloc( N*sizeof(double*) );
  for( i=0; i<N; i++){
	 d[i] = (double*) malloc( M*sizeof(double) );
	 for( j=0; j<M; j++ ){
		d[i][j]=0.0;
	 }
  }
  
  return d;
}
/** return a 0-initialized integer matrix of dimension NxM
 */
int** matrix_init_int(int N, int M){
  int i,j;
  int **d;
  /* dprintf("N,M=(%i,%i)\n", N, M); */
  d = (int**) malloc( N*sizeof(int*) );
  for( i=0; i<N; i++){
 	 /* dprintf(" i=%i\n", i ); */
	 d[i] = (int*) malloc( M*sizeof(int) );
	 for( j=0; j<M; j++ ){
		d[i][j]=0;
	 }
  }
  
  return d;
}


/** Delete a column in a matrix. Memory remains allocated and everything is 
 *  moved (deleted column sits at index n). Index runs from 0,...,n-1
 * \param m - matrix 
 * \param N - num of rows
 * \param n - num of cols
 * \param col- num of col to delete
 */
double** matrix_delcol(double **m, int N, int n, int col){
  int i, j;
  double tmp;
  massert(col>n-1, "cannot delete col %i because only %i cols\n", col, n);
  for(i=0; i<N; i++){
    tmp = m[i][col];
    for(j=col; j<n-1; j++){
      m[i][j]=m[i][j+1];
    }
    m[i][n-1]=tmp;
  }
  return m;
}

/** Get minimum entry from matrix. Set indices accordingly.
 * \param m, N, n - Nxn matrix
 * \param i1,i2 - pointer to indices of min-element. If NULL, ignored.
 */
double matrix_min(const double **m, int N, int n, int *i1, int *i2){
  int i,j;
  double minel=DBL_MAX;
  for(i=0; i<N; i++){
    for(j=0; j<n; j++){
      if(m[i][j]<minel){
		  minel = m[i][j];
		  if(i1 && i2){
			 *i1 = i;
			 *i2 = j;
		  }
      }
    }
  } 

  if(i1 && i2){
	 dprintf("matrix_min: m(%i,%i)=%f (%f)\n", *i1, *i2, m[*i1][*i2], minel);
  } else {
	 dprintf("matrix_min: m=%f\n",  minel);
  }
  return minel;
}
/** Get maximum entry from matrix. Set indices accordingly.
 * \param m, N, n - Nxn matrix
 * \param i1,i2 - pointer to indices of max-element. If NULL, ignored.
 */
double matrix_max(const double **m, int N, int n, int *i1, int *i2){
  int i,j;
  double maxel=DBL_MIN;
  for(i=0; i<N; i++){
    for(j=0; j<n; j++){
      if(m[i][j]>maxel){
		  maxel = m[i][j];
		  if(i1 && i2){
			 *i1 = i;
			 *i2 = j;
		  }
      }
    }
  }
  if(i1 && i2){
	 dprintf("matrix_max: m(%i,%i)=%f (%f)\n", *i1, *i2, m[*i1][*i2], maxel);
  } else {
	 dprintf("matrix_max: %f\n", maxel);
  }
  return maxel;
}


/** print a matrix.
 */
void matrix_print(double **m, int N, int n){
  int i,j;

  for(i=0; i<N; i++){
	 printf("[ ");
	 for(j=0; j<n; j++){
		if(m[i][j]>1000000)
		  printf("<max> ");
		else
		  printf("%2.2f ", m[i][j]);
	 }
	 printf(" ]\n");
  }
}
/** add matrix m2 to matrix m1 (must be of equal dimensions
 */
void matrix_add_matrix(double **m1, const double **m2, int N, int n){
  int i, j;

  for( i=0; i<N; i++){
	 for( j=0; j<n; j++ ){
		m1[i][j]+=m2[i][j];
	 }
  }
}

/** 
	 MatLab's M1.*M2, which is M1[i][j] = M1[i][j]*M2[i][j];
	 That is, m1 is overwritten with the result from the operation
	 \param m1,m2
	 \param N,M dimensions of m1 and m2 
 */
void     matrix_dottimes_matrix( double **m1, const double **m2, int N, int M ){
  int i,j;

  for( i=0; i<N; i++){
	 for( j=0; j<M; j++ ){
		m1[i][j] = m1[i][j]*m2[i][j];
	 }
  }
}

/** copy matrix src to dest with dimensions N,M
 */
void     matrix_copy( const double **src, double **dest, int N, int M ){
  int i,j;

  for( i=0; i<N; i++ ){
	 for( j=0; j<M; j++ ){
		dest[i][j] = src[i][j];
	 }
  }
}
/** divide all entries in m by s
 */
void matrix_divide_scalar(double **m, int N, int n, double s){
  int i, j;

  for( i=0; i<N; i++){
	 for( j=0; j<n; j++ ){
		m[i][j] /= s;
	 }
  }
}
void     matrix_add_scalar(double **m, int N, int n, double s){
  int i, j;

  for( i=0; i<N; i++){
	 for( j=0; j<n; j++ ){
		m[i][j] += s;
	 }
  }
}
/** multiply all entries in m by s
 */
void matrix_mul_scalar(double **m, int N, int n, double s){
  int i, j;

  for( i=0; i<N; i++){
	 for( j=0; j<n; j++ ){
		m[i][j] *= s;
	 }
  }
}

/** free matrix memory 
 */
void matrix_free(double **m, int N){
  int i;
  if( !m ){
	 return;
  }
  for( i=0; i<N; i++){
	 free(m[i]);
  }
  free(m);
}

/** signum function 
 */
int sgn(int x){
  return (x > 0) ? 1 : (x < 0) ? -1 : 0;
}


/** Compute how many points are in the discretized straight line computed
	 by Bresenham's line-drawing algorithm. Can be used to determine how much
	 memory needs to be allocated.
 *  \param  int xstart, ystart        = Koordinaten des Startpunkts
 *  \param  int xend, yend            = Koordinaten des Endpunkts
 */
int bresenham_howmany_points( int xstart,int ystart,int xend,int yend ){
  int dx, dy;
  dx = (xend - xstart);
  dy = (yend - ystart);
  return MAX( dx, dy );
}

/**
 * Bresenham-Algorithm (stolen but modified from Wikipedia)
 *
 *  \param  int xstart, ystart        = Koordinaten des Startpunkts
 *  \param  int xend, yend            = Koordinaten des Endpunkts
 *  \param points -- output; saved in pairs (x,y) of the coordinates;
                     must contain enough allocated memory 2*(maximum of 
							(yend-ystart) and (xend-xstart); if NULL is passed, 
							memory is alloced within the function
*/
int* bresenham(int xstart,int ystart,int xend,int yend, int *points){
  int x, y, t, dx, dy, incx, incy, pdx, pdy, ddx, ddy, es, el, err;
  int points_idx;
  points_idx = 0;

  /* Entfernung in beiden Dimensionen berechnen */
  dx = xend - xstart;
  dy = yend - ystart;
 
  if( points==ALLOC_IN_FCT ){
	 warnprintf( "allocating in bresenham!\n");
	 points = (int*) malloc( 2*bresenham_howmany_points( xstart, ystart, 
																		  xend, yend )*sizeof(int) );
  }

  /* Vorzeichen des Inkrements bestimmen */
  incx = sgn(dx);
  incy = sgn(dy);
  if(dx<0) dx = -dx;
  if(dy<0) dy = -dy;
 
  /* feststellen, welche Entfernung größer ist */
  if (dx>dy)
	 {
      /* x ist schnelle Richtung */
      pdx=incx; pdy=0;    /* pd. ist Parallelschritt */
      ddx=incx; ddy=incy; /* dd. ist Diagonalschritt */
      es =dy;   el =dx;   /* Fehlerschritte schnell, langsam */
	 } else
	 {
      /* y ist schnelle Richtung */
      pdx=0;    pdy=incy; /* pd. ist Parallelschritt */
      ddx=incx; ddy=incy; /* dd. ist Diagonalschritt */
      es =dx;   el =dy;   /* Fehlerschritte schnell, langsam */
	 }
 
  /* Initialisierungen vor Schleifenbeginn */
  x = xstart;
  y = ystart;
  err = el/2;
  points[0]=x;
  points[1]=y;	
 
  /* Pixel berechnen */
  for(t=1; t<=el; t++) /* t zaehlt die Pixel, el ist auch Anzahl */
	 {
      /* Aktualisierung Fehlerterm */
      err -= es; 
      if(err<0)
		  {
          /* Fehlerterm wieder positiv (>=0) machen */
          err += el;
          /* Schritt in langsame Richtung, Diagonalschritt */
          x += ddx;
          y += ddy;
		  } else
		  {
          /* Schritt in schnelle Richtung, Parallelschritt */
          x += pdx;
          y += pdy;
		  }
		points[( 2*t )+0]=x;
		points[( 2*t )+1]=y;	
		/* dprintf(" 2*t+1=%i\n", 2*t+1); */
	 }
  //	dprintf("t=%i, (x,y)=(%i,%i)\n", t,x,y);
  return points;
} /* bresenham() */

void    swap2i(int *v1, int *v2){
  int tmp;
  tmp = *v1;
  *v1 = *v2; 
  *v2 = tmp;
}

void    swap2d(double *v1, double *v2){
  double tmp;
  tmp = *v1;
  *v1 = *v2; 
  *v2 = tmp;
}



/** Draw a sample from a discrete distribution (nearest neighbour).
	 \param v,n the distribution as a funciton of [0,...,n]
	 \param x a real value from the interval [0,...,n]
	 \return v(x) by nn-interpolation
*/
double drawsample_nearest_neighbour( const double *v, int n, double x ){
  int x1;
  double r;

  x1 = (int)round( x );

  if( x1<0 || x1>=n ){
	 errprintf( "(x=%f, x1=%i, n=%i)\n", x, x1, n );
  }
  r = v[x1];
  if( isnan( r ) ){
	 errprintf( "r=%f (x=%f, x1=%i,v[x1]=%f, n=%i)\n", r, x, x1,v[x1],n );
  }
  return r;
}


/** Resample a vector to match a new length (NN-interpolation). 
	 Assume that s is from [0,...,n] and resample, such that the new
	 x is running from [0,...,newn] where y[newn]=y_old[n].
	 \param s,n the original vector
	 \param newn the new length of the vector
	 \param news caller-allocated memory of at least length newn (if NULL, 
	             the function allocates own memory)
	 \return resampled vector in news
	 \see drawsample_nearest_neighbour()
*/
double* resample_nearest_neighbour( const double *s, int n, int newn, double *news ){
  double step;
  int i;

  if( news==NULL ){
	 news = (double*) malloc( newn*sizeof( double ) );
  }

  step = (double)n/(double)newn;
  for( i=0; i<newn; i++ ){
	 news[i] = drawsample_nearest_neighbour( s, n, i*step );
  }
 
  return news;
}

/** Resample a vector to match a new length (linear interpolation).
	 \param s,n the original vector
	 \param newn the new length of the vector
	 \param news caller-allocated memory of at least length newn (if NULL, 
	             the function allocates own memory)
	 \return resampled vector in news
	 \see drawsample_linear()
*/
double* resample_linear( const double *s, int n, int newn, double *news ){
  return resample_gsl( s, n, newn, news, gsl_interp_linear );
}

/** Resample a vector to match a new length using GSL's interpolation
	 options:
	  -# gsl_interp_linear
	  -# gsl_interp_polynomial
	  -# gsl_interp_cspline
	  -# gsl_interp_cspline_periodic
	  -# gsl_interp_akima
	  -# gsl_interp_akima_periodic

	 Assume that s is from [0,...,n] and resample, such that the new
	 x is running from [0,...,newn] where y[newn]=y_old[n].
	 \param s,n the original vector
	 \param newn the new length of the vector
	 \param news caller-allocated memory of at least length newn (if NULL, 
	             the function allocates own memory)
    \param method one of GSL interpolation types (see above)
	 \return resampled vector in news
*/
double* resample_gsl( const double *s, int n, int newn, double *news, const gsl_interp_type *method ){
  double step;
  int i;
  double *x;

  x = (double*)malloc( n*sizeof(double) );
  for( i=0; i<n; i++ ){
	 x[i] = (double)i;
  }
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (method, n);
  gsl_spline_init (spline, x, s, n);

  if( news==NULL ){
	 news = (double*) malloc( newn*sizeof( double ) );
  }


  step = (double)n/(double)newn;
  for( i=0; i<newn; i++ ){
	 news[i] = gsl_spline_eval( spline, i*step, acc );
  }

  /* free */
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  free(x);

  return news;
}
double* flip_array( double *v, int n ){
  int i;
  double tmp;

  for( i=0; i<n/2; i++){
	 tmp = v[i];
	 v[i] = v[n-i-1];
	 v[n-i-1] = tmp;
  }
  return v;
}



#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
/** 1D Fourier transform (FFT). Taken from Numerical Recipes in C.
	 \param data - array to be tranformed (input/output); for input, 
	               is assumed to be complex, i.e. [ real, imag, real, imag, ... ]
						data contains nn complex numbers and 2*nn entries.
	 \param nn - must be a power of two, not checked for!; length(data)=2*nn
	 \param isign - 1: forward; -1: ifft;
	 \return data is now the FT in complex numbers
 */
void fft(double *data, unsigned long nn, int isign){
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  float tempr,tempi;
  data--; /* because num_rec assumes [1,...n] arrays */
  /* dprintf( "nn=%li, isign=%i\n", nn, isign ); */

  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
	 if (j > i) {
		SWAP(data[j],data[i]);
		SWAP(data[j+1],data[i+1]);
	 }
	 m=nn;
	 while (m >= 2 && j > m) {
		j -= m;
		m >>= 1;
	 }
	 j += m;
  }
  mmax=2;
  while (n > mmax) {
	 istep=mmax << 1;
	 theta=isign*(6.28318530717959/mmax);
	 wtemp=sin(0.5*theta);
	 wpr = -2.0*wtemp*wtemp;
	 wpi=sin(theta);
	 wr=1.0;
	 wi=0.0;
	 for (m=1;m<mmax;m+=2) {
		for (i=m;i<=n;i+=istep) {
		  j=i+mmax;
		  tempr=wr*data[j]-wi*data[j+1];
		  tempi=wr*data[j+1]+wi*data[j];
		  data[j]=data[i]-tempr;
		  data[j+1]=data[i+1]-tempi;
		  data[i] += tempr;
		  data[i+1] += tempi;
		}
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	 }
	 mmax=istep;
  }
}
#undef SWAP

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
/** \code [1 2 3 - - - - -] -> [1 1 1 2 3 3 3 3] \endcode */
double* sigext_smooth(double *data, int ns, int n){
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

/** return min_k( 2^k | 2^k>n ).
 */ 
int next_pow2( int n ){
  int p;
  p = (int)round(glog((double)n, 2))+1;
  return (int)p;//(int)pow(2, p);
}

/** computes remainder of x/y using
	 \f$ 
	 r(x,y) = \left\lfloor x - y\left\lfloor \frac{x}{y}\right\rfloor \right\rfloor
	 \f$
 */
int iremainder( double x, double y){
  int result;
  
  if (y != 0) {
	 result =  x-y* (int)(x/y);
  } else  {
	 result = 0;
	 errprintf("division by zero!\n");
  }

  return result;
}

double** disttransform_deadreckoning(int **I, int X, int Y, double **d){
   int x, y;
   double d1=1;
   double d2=sqrt(2);
   int **P1, **P2;

	P1 = matrix_init_int( X, Y );
	P2 = matrix_init_int( X, Y );

   /* P1 = (int**) malloc(X*sizeof(int*)); */
   /* P2 = (int**) malloc(X*sizeof(int*)); */
	/* printf("fdas, %i, %i\n", X, Y); */
   /* for(x=0; x<X; x++){ */
	/*   printf("1, %i, %i\n", x, Y); */
	/*   P1[x]=(int*)calloc(Y, sizeof(int)); */
	/*   printf("2, %i, %i\n", x, Y); */
	/*   P2[x]=(int*)calloc(Y, sizeof(int)); */
   /* } */
	/* printf("fdas\n"); */

	if( d==NULL ){
	  d = matrix_init( X, Y );
	}
  
   for(y=0; y<Y; y++){
	  for(x=0; x<X; x++){
		 d[x][y] = DBL_MAX;
		 P1[x][y]=-1;
		 P2[x][y]=-1;
		 if(x>0 && y>0 && y<Y-1 && x<X-1){
			if( I[x-1][y  ]!=I[x][y] || I[x+1][y  ]!=I[x][y] ||
				 I[x  ][y-1]!=I[x][y] || I[x  ][y+1]!=I[x][y]){
			  d [x][y]=0;
			  P1[x][y]=x;
			  P2[x][y]=y;
			}
		 }            
	  }
   }
   
   /* first pass */
   for(y=1; y<Y-1; y++){
	  for(x=1; x<X-1; x++){
		 if(d[x-1][y-1]+d2 < d[x][y]){
			P1[x][y] = P1[x-1][y-1];
			P2[x][y] = P2[x-1][y-1];
			d [x][y] = sqrt( SQR( x-P1[x][y] ) + 
								  SQR( y-P2[x][y] ) );
		 }
		 if(d[x  ][y-1]+d1 < d[x][y]) {
			P1[x][y] = P1[x][y-1];
			P2[x][y] = P2[x][y-1];
			d[x][y]  = sqrt( SQR( x-P1[x][y] ) + 
								  SQR( y-P2[x][y] ) );
                
		 }           
		 if(d[x+1][y-1]+d2 < d[x][y]) {
			P1[x][y] = P1[x+1][y-1];
			P2[x][y] = P2[x+1][y-1];
			d[x][y]  = sqrt( SQR( x-P1[x][y] ) + 
								  SQR( y-P2[x][y] ) );
                
		 }
		 if(d[x-1][y  ]+d1 < d[x][y]) {
			P1[x][y] = P1[x-1][y];
			P2[x][y] = P2[x-1][y];
			d[x][y]  = sqrt( SQR( x-P1[x][y] ) + 
								  SQR( y-P2[x][y] ) );
                
		 }
	  }
   }

   /* final pass */
   for(y=Y-2; y>0; y--){
	  for(x=X-2; x>0; x--){
		 if(d[x+1][y]+d1 < d[x][y]){
			P1[x][y] = P1[x+1][y];
			P2[x][y] = P2[x+1][y];
			d[x][y]  = sqrt( SQR( x-P1[x][y] ) + 
								  SQR( y-P2[x][y] ) );
		 }
		 if(d[x-1][y+1]+d2 < d[x][y]) {
			P1[x][y] = P1[x-1][y+1];
			P2[x][y] = P2[x-1][y+1];
			d[x][y]  = sqrt( SQR( x-P1[x][y] ) + 
								  SQR( y-P2[x][y] ) );
                
		 }           
		 if(d[x][y+1]+d1 < d[x][y]) {
			P1[x][y] = P1[x][y+1];
			P2[x][y] = P2[x][y+1];
			d[x][y]  = sqrt( SQR( x-P1[x][y] ) + 
								  SQR( y-P2[x][y] ) );
                
		 }
		 if(d[x+1][y+1]+d2 < d[x][y]) {
			P1[x][y] = P1[x+1][y+1];
			P2[x][y] = P2[x+1][y+1];
			d[x][y]  = sqrt( SQR( x-P1[x][y] ) + 
								  SQR( y-P2[x][y] ) );
                
		 }
	  }
   }
    
   for(x=0; x<X; x++){
	  d[x][0  ] = d[x][1];
	  d[x][X-1] = d[x][X-2];
   }
   for(y=0; y<Y; y++){
	  d[0  ][y] = d[1  ][y];
	  d[Y-1][y] = d[Y-2][y];
   }
   for(y=Y-2; y>0; y--)
	  for(x=X-2; x>0; x--)
		 if(I[x][y]>0) d[x][y] = -1*d[x][y];
   

   for(x=0; x<X; x++){
	  free(P1[x]);
	  free(P2[x]);
   }
   free(P1); free(P2);
   
	return d;
}

/** \f$
	 G(x) = \frac{1}{\sigma\sqrt{2\pi}} \exp{\left(-\frac{(x-\mu)^2}{2\sigma^2}\right)}
	 \f$
*/
double gaussian( double x, double sigma, double mu ){
  return  1/(sigma*sqrt(2*PI)) * exp( -1 * ( SQR( (x-mu) ) / 
															(2*SQR( sigma )) ) );
}


/** Computes: m[i][j] = scalar-m[i][j]
 */
void scalar_minus_matrix( double scalar, double **m, int N, int M ){
  int i,j;


  for( i=0; i<N; i++ ){
	 for( j=0; j<M; j++ ){
		m[i][j] = scalar-m[i][j];
	 }
  }
}


Complex complex( double re, double im ){
  Complex a;
  a.re = re;
  a.im = im;
  return a;
}
Complex complex_add( Complex a, Complex b ){
  a.re += b.re;
  a.im += b.im;
  return a;
}
Complex complex_sub( Complex a, Complex b ){
  a.re -= b.re;
  a.im -= b.im;
  return a;
}
Complex complex_mul( Complex a, Complex b ){
  Complex r;
  r.re = (a.re*b.re)-(a.im*b.im);
  r.im = (a.im*b.re)+(a.re*b.im);
  return r;
}
Complex complex_mul_double( Complex a, double b ){
  a.re *= b;
  a.im *= b;
  return a;
}
double  complex_abs( Complex a ){
  return sqrt( SQR( a.re ) + SQR( a.im ) );
}

/** this is
	 \f[
	 \exp( a+bi ) = \exp( a )\cdot (\cos b + i\sin b )
	 \f]
 */
Complex complex_exp( Complex a ){
  Complex r;
  r = complex( cos( a.im ), sin( a.im ) );
  r = complex_mul_double( r, exp( a.re ) );
  return r;
}
Complex complex_conj( Complex a ){
  a.im = -a.im;
  return a;
}

Complex complex_neg ( Complex a ){
  a.re = -a.re;
  a.im = -a.im;
  return a;
}

Complex complex_sqrt( Complex x ){
  double r = complex_abs(x);
  Complex z = complex(sqrt(0.5 * (r + x.re)),
							 sqrt(0.5 * (r - x.re)));
  if (x.im < 0.0) z.im = -z.im;
  return z;
}
/** this is
\f[ 
\frac{(a + bi)}{(c + di)} = \left({ac + bd \over c^2 + d^2}\right) + \left( {bc - ad \over c^2 + d^2} \right)i
\f]
*/
Complex complex_div( Complex a, Complex b){
  Complex r;
  r.re = (a.re*b.re + a.im*b.im)/( SQR( b.re )+SQR( b.im ) );
  r.im = (a.im*b.re - a.re*b.im)/( SQR( b.re )+SQR( b.im ) );
  return r;
}

Complex complex_bilinear_transform(Complex pz){
  return complex_div( complex_add( complex(2.0,0.0), pz ), complex_sub( complex(2.0,0.0), pz ) );
}

/** expand a polynomial given in root-form
	 \f[
	 p(x) = x_0(x-x_1)(x-x_2)\cdots (x-x_n)
	 \f]
	 to find the coefficients a_0, ..., a_n of
	 \f[
	 p(x) = a_0 + a_1 x_1 + \cdots + a_n x_n
	 \f]
	\param roots x_0,...,x_n
	\param coeffs n+1 numbers to hold the resulting coefficients or ALLOC_IN_FCT
*/
Complex* expand_polynomial_from_roots( const Complex *roots, int n, Complex *coeffs ){
  int i, j;
  Complex nw;

  if( coeffs==ALLOC_IN_FCT ){
	 coeffs = (Complex*) malloc( (n+1)*sizeof( Complex ) );
  }

  coeffs[0] = complex( 1.0, 0.0 );
  for (i=0; i < n; i++){
	 coeffs[i+1] = complex(0.0,0.0);
  }
  for (i=0; i < n; i++){
	 nw = complex_mul_double( roots[i], -1.0 );
	 for( j=n; j>=1; j-- ){
		coeffs[j] = complex_mul( nw, coeffs[j] );
		coeffs[j] = complex_add( coeffs[j], coeffs[j-1] );
	 }
	 coeffs[0] = complex_mul( nw, coeffs[0] );
  }


  return coeffs;
}


/** return the real part of the complex numbers in vc.
	 If the imaginary part is nonzero, a warning is issued.
	 \param vc complex vector
	 \param vr real vector of size n, or ALLOC_IN_FCT
*/
double* vector_complex_to_real( const Complex *vc, double *vr, int n ){
  int i;
  dprintf("entering\n");
  if( vr==ALLOC_IN_FCT ){
	 warnprintf("allocating in fct\n");
	 vr = vector_init( vr, n, 0.0 );
  }
  for( i=0; i<n; i++ ){
	 vr[i] = vc[i].re;
	 if( ABS(vc[i].im)>1e-10 ){
		warnprintf("complex number with index '%i' has nonzero imaginary part\n", i);
	 }
  }
  return vr;
}
