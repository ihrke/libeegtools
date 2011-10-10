/*
** check_list.c
** 
** Unit-Test suite using libcheck
**
**
** Made by (Matthias Ihrke)
** Login   <mihrke@localhost>
** 
*/

#include<stdlib.h>
#include<check.h>
#include<math.h>
#include "check_all.h"
#include "helper.h"
#include "optarg.h"
#include "eeg.h"
#include "mpitools.h"
#include "pqueue.h"
#include "mathadd.h"
#include "nnsearch.h"
#include "imageproc.h"
#include "time_frequency.h"
#include "io.h"
#include "array.h"
#include "linalg.h"
     
START_TEST (test_strip_blank)
{
  char s1[40]=" thi s is a blank st ri n g ";
  char s2[40]="  more blanks     at the end     ";

  string_strip_blanks( (char*)s1 );
  string_strip_blanks( (char*)s2 );

  fail_if( strcmp( s1, "thisisablankstring" ), "FAILURE: got '%s'\n", s1 );
  fail_if( strcmp( s2, "moreblanksattheend" ), "FAILURE: got '%s'\n", s2 );

}
END_TEST


START_TEST (test_eeg_to_stream)
{

}
END_TEST

START_TEST (test_nnprepare)
{
#ifdef EXPERIMENTAL
  double **X;
  double *q;
  int N, p;
  int i,j;
  N=51;
  p=3;
  X = (double**)malloc( N*sizeof(double*) );
  q = (double*) malloc( p*sizeof(double ) );
  for( i=0; i<N; i++){
	 X[i] = (double*)malloc( p*sizeof(double*) );
	 for( j=0; j<p; j++ ){
		X[i][j] = drand48();
	 }
  }
  for( i=0; i<p; i++ ){
	 q[i] = drand48();
  }

  OptArgList *o;
  o=optarglist("max_numel_terminal_node=int", 2);
  optarglist_print( o, stderr );

  SearchTree *S;
  S = nn_prepare( X, N, p, o );

  dprintf("S->N=%i\n", S->N );
  int k=N/5;
  int* nn_idx = (int*)malloc( k*sizeof(int));
  double *nnd = (double*)malloc( k*sizeof(double));
  dprintf("looking for k=%i neighbours of vector\n",k);
  dblp_print( q, p );

  nn_search_k( S, q, k, nn_idx, nnd );
  dprintf("Finished searching, result: \n");
  dblp_print( nnd, k );
  dblp_print_int(nn_idx, k);

  dprintf("slow searching:\n");
  //  nn_search_k_slow( X, N, p, q, k, nn_idx, nnd, NULL );
  dprintf("Finished searching, result: \n");
  dblp_print( nnd, k );
  dblp_print_int(nn_idx, k);

  for( i=0; i<N; i++){
	 free( X[i] );
  }
  free( X );
  free( q );
  free( nn_idx );
  free( nnd );
#endif					
}
END_TEST

START_TEST (test_pqueue)
{
  PriorityQueue *pq=pq_init();
  int stuff[] = {1, 2, 3, 4};
  int i;
  dprintf("Inserting into queue\n");
  pq_insert( pq, (void*)stuff, 10 );
  pq_insert( pq, (void*)(stuff+1), 5 );
  pq_insert( pq, (void*)(stuff+2), 7 );
  pq_insert( pq, (void*)(stuff+3), 1 );
  dprintf("Done\n");
  pq_print( stderr, pq );
  int p[4];
  for( i=0; i<4; i++ ){
	 p[i]=*((int*)pq_pop( pq ));
	 dprintf("popped %i\n", p[i]);
  }
  fail_unless( stuff[0]==p[0] );
  fail_unless( stuff[1]==p[2] );
  fail_unless( stuff[2]==p[1] );
  fail_unless( stuff[3]==p[3] );

  pq_free( pq );
}
END_TEST

START_TEST (test_bresenham)
{
  int as[2]={2,3};
  int ae[2]={10,12};
  Array *a=bresenham_line( as, ae );
  fail_unless( a->ndim==2 );
  fail_unless( a->size[0]==2 );
  fail_unless( a->dtype==INT );
  /* start */
  fail_unless( array_INDEX2(a,int,0,0)==as[0] &&
					array_INDEX2(a,int,1,0)==as[1] );
  /* end */
  fail_unless( array_INDEX2(a,int,0,a->size[1]-1)==ae[0] &&
					array_INDEX2(a,int,1,a->size[1]-1)==ae[1] );
					
  int i; /* coninuity */
  for( i=1; i<a->size[1]; i++ ){
	 fail_if( array_INDEX2(a,int,0,i)-array_INDEX2(a,int,0,i-1)>1 ||
				 array_INDEX2(a,int,0,i)-array_INDEX2(a,int,0,i-1)<0 );
	 fail_if( array_INDEX2(a,int,1,i)-array_INDEX2(a,int,1,i-1)>1 ||
				 array_INDEX2(a,int,1,i)-array_INDEX2(a,int,1,i-1)<0 );
  }
  array_free( a );
  
  int bs[2]={-100,2};
  int be[2]={200,30};
  Array *b=bresenham_line( bs, be );
  
  fail_unless( b->ndim==2 );
  fail_unless( b->size[0]==2 );
  fail_unless( b->dtype==INT );


  /* start */
  fail_unless( array_INDEX2(b,int,0,0)==bs[0] &&
					array_INDEX2(b,int,1,0)==bs[1] );
  /* end */
  fail_unless( array_INDEX2(b,int,0,b->size[1]-1)==be[0] &&
					array_INDEX2(b,int,1,b->size[1]-1)==be[1] );

  /* coninuity */
  for( i=1; i<b->size[1]; i++ ){
	 fail_if( array_INDEX2(b,int,0,i)-array_INDEX2(b,int,0,i-1)>1 ||
				 array_INDEX2(b,int,0,i)-array_INDEX2(b,int,0,i-1)<0 );
	 fail_if( array_INDEX2(b,int,1,i)-array_INDEX2(b,int,1,i-1)>1 ||
				 array_INDEX2(b,int,1,i)-array_INDEX2(b,int,1,i-1)<0 );
  }

  array_free( b );
}
END_TEST

START_TEST (test_bresenham_segments)
{
  /*
	 load bres.mat
	 x=zeros(200,200);
	 for i=1:size(bres,2)
	    x(bres(1,i),bres(2,i))=1;
	 end;
	 imagesc(x);
	 colormap gray;
	*/
  int points[12]={ 3, 7, 100, 50, 10, 150, 
						 4, 10, 20, 100, 110, 200 };
  Array *p=array_fromptr2( INT, 2, points, 2, 6 );
  Array *a=bresenham_linesegments( p );

  fail_unless( a->ndim==2 );
  fail_unless( a->size[0]==2 );
  fail_unless( a->dtype==INT );

  int i; /* coninuity */
  for( i=1; i<a->size[1]; i++ ){
	 fail_if( array_INDEX2(a,int,0,i)-array_INDEX2(a,int,0,i-1)>1 ||
				 array_INDEX2(a,int,0,i)-array_INDEX2(a,int,0,i-1)<0 );
	 fail_if( array_INDEX2(a,int,1,i)-array_INDEX2(a,int,1,i-1)>1 ||
				 array_INDEX2(a,int,1,i)-array_INDEX2(a,int,1,i-1)<0 );
  }


  /* write_array_matlab( a, "bres", "bres.mat", FALSE );  */

  array_free( a );
  array_free( p );
}
END_TEST

START_TEST (test_disttransform)
{
  /*
	 load test.mat
	 subplot(1,2,1)
	 imagesc(in)
	 colorbar
	 subplot(1,2,2)
	 imagesc(dt)
	 colorbar
	 colormap gray
  */
  int nx=100, ny=200;
  Array *a=array_new2(INT,2,nx,ny);
  array_INDEX2(a,int,nx/2,ny/2)=1;
  array_INDEX2(a,int,1,1)=1;
  array_INDEX2(a,int,1,ny-1)=1;

  Array *b=disttransform_deadreckoning( a, NULL );

  ulong i; 
  for( i=0; i<array_NUMEL(b); i++ ){
	 fail_if( isnan( array_INDEX1(b,double,i) ) );
	 fail_if( array_INDEX1(b,double,i)>DBL_MAX-2 );
  }

  int as[2]={2,3};
  int ae[2]={100,120};
  Array *l=bresenham_line( as, ae );
  for( i=0; i<l->size[1]; i++ ){
	 array_INDEX2(a,int, array_INDEX2(l,int,0,i),
					  array_INDEX2(l,int,1,i) )=1;
  }
  Array *d=disttransform_deadreckoning( a, NULL );
 for( i=0; i<array_NUMEL(b); i++ ){
	 fail_if( isnan( array_INDEX1(d,double,i) ) );
	 fail_if( array_INDEX1(d,double,i)>DBL_MAX-2 );
  }
  /* write_array_matlab( a, "in", "test.mat", FALSE ); */
  /* write_array_matlab( d, "dt", "test.mat", TRUE ); */

  array_free( a );
  array_free( b );
  array_free( l );
  array_free( d );  
}
END_TEST

START_TEST (test_regularization_line)
{
  /*
	 load reg.mat
	 imagesc(reg);
	*/
  int points[12]={ 0, 50, 60, 65, 150, 199, 
						 0, 10, 20, 100, 110, 199 };
  int dims[2]={200,200};
  Array *p=array_fromptr2( INT, 2, points, 2, 6 );
  Array *a=regularization_linear_points( p, dims, NULL );

  bool ismat;
  matrix_CHECK(ismat,a);
  fail_unless( ismat );

  /* write_array_matlab( a, "reg", "reg.mat", FALSE );   */

  array_free( a );
  array_free( p );
}
END_TEST
START_TEST (test_regularization_gauss)
{
  /*
	 load reg.mat
	 imagesc(reg);
	*/
  int points[12]={ 0, 50, 60, 65, 150, 199, 
						 0, 10, 20, 100, 110, 199 };
  int dims[2]={200,200};
  Array *p=array_fromptr2( INT, 2, points, 2, 6 );
  Array *a=regularization_gaussian( p, dims, NULL, 0.2 );

  bool ismat;
  matrix_CHECK(ismat,a);
  fail_unless( ismat );

  /* write_array_matlab( a, "reg", "reg.mat", FALSE );*/

  array_free( a );
  array_free( p );
}
END_TEST
START_TEST (test_regularization_gauss_narrowdown)
{
	/*
	 load reg.mat
	 imagesc(reg);
	*/
  int points[12]={ 0, 50, 60, 65, 150, 199,
				 0, 10, 20, 100, 110, 199 };
  int dims[2]={200,200};
  Array *p=array_fromptr2( INT, 2, points, 2, 6 );
  Array *a=regularization_gaussian_narrowdown( p, dims, NULL, 0.4 );

  bool ismat;
  matrix_CHECK(ismat,a);
  fail_unless( ismat );

  /* write_array_matlab( a, "reg", "reg.mat", FALSE );*/

  array_free( a );
  array_free( p );
}
END_TEST

START_TEST (test_spectrogram_init)
{
  int i,j;
  Spectrogram *s=spectrogram_init( 100, 1000 );
  fail_unless( s->N_freq==100 );
  fail_unless( s->N_time==1000 );
  for( i=0; i<1000; i++ ){
	 for( j=0; j<100; j++ ){
		s->sgram[i][j].re;
		s->sgram[i][j].im;
	 }
  }
  spectrogram_free( s );
}
END_TEST

START_TEST (test_windows)
{
  int n=100;
  bool isvec;

  Array *w1=window_dirichlet( n, 0.0 );
  vector_CHECK( isvec,w1);
  fail_unless( isvec );
  Array *w2=window_gaussian( n, 0.3 );
  vector_CHECK( isvec,w2);
  fail_unless( isvec );
  Array *w3=window_hamming( n, 0.0 );
  vector_CHECK( isvec,w3);
  fail_unless( isvec );
  Array *w4=window_hanning( n, 0.0 );
  vector_CHECK( isvec,w4);
  fail_unless( isvec );
  Array *w5=window_kaiser( n, 1.0 );
  vector_CHECK( isvec,w5);
  fail_unless( isvec );

  array_free( w1 );
  array_free( w2 );
  array_free( w3 );
  array_free( w4 );
  array_free( w5 );

}
END_TEST

START_TEST (test_spectgram)
{
  /*
	 load spec.mat
	 [Nt,Nf]=size(powspec);
	 fx=linspace(0,500,Nf);
	 tx=linspace(0,1,Nt);
	 imagesc(tx,fx,powspec);
	 ylabel('Frequency (Hz)');
	 xlabel('Time (arbitrary)');
	*/
#ifdef MATIO
  Array *a=read_array_matlab( CHECKDATADIR"/chirp.mat", NULL );
  array_dimred( a ); /* MATLAB writes as 1 x N */
  /* array_print( a, -1, stderr );  */

  bool isvec;
  vector_CHECK(isvec,a);
  fail_if( !isvec);
  
  Array *w=window_hamming( 129, 0.0 );
  double srate=1000.0;
  int Nf=128;
  int Nt=a->size[0]/10;
  double cf[2]={0.0, 500.0};
  Spectrogram *s= spectrogram_stft(a, srate,w,Nf, Nt,cf, NULL, NULL );
  Array *p=spectrogram_powerspectrum( s );

  bool ismat;
  matrix_CHECK(ismat,p);
  fail_if( !ismat);
  fail_if( p->size[0]!=s->N_time );
  fail_if( p->size[1]!=s->N_freq );
  /* write_array_matlab( p, "powspec", "spec.mat", FALSE );  */

  spectrogram_free( s );
  array_free( a );
  array_free( p );
  array_free( w );
#endif
}
END_TEST

START_TEST (test_check_bit)
{
  fail_unless( CHECK_BIT( 1, 0 ) );
  fail_unless( !CHECK_BIT( 1, 1 ) );

  fail_unless( CHECK_BIT( 7, 0 ) );
  fail_unless( CHECK_BIT( 7, 1 ) );
  fail_unless( CHECK_BIT( 7, 2 ) );

  fail_unless( !CHECK_BIT( 4, 0 ) );
  fail_unless( !CHECK_BIT( 4, 1 ) );
  fail_unless( CHECK_BIT( 4, 2 ) );
}
END_TEST
START_TEST (test_interp1)
{
	/* MATLAB
	   load interp1.mat
	   plot( x, y, 'rx-' );
	   hold on
	   plot( xi, yi, 'b-o' );
	   hold off
	   legend('original', 'interpolated');
	   */
	Array *x=array_new_dummy( DOUBLE, 1, 10 );
	Array *y=array_new2( DOUBLE, 1, 10 );
	int i;
	for( i=0; i<10; i++ ){
		vec_IDX( y, i)=sin(vec_IDX(x,i));
	}
	Array *xi=array_new2( DOUBLE, 1, 100 );
	for( i=0; i<100; i++ ){
		vec_IDX(xi, i)=((double)i)/100.0*10.0;
	}
	Array *yi=vector_interp1( x, y, xi, NULL );
	write_array_matlab( x, "x", "interp1.mat", FALSE );
	write_array_matlab( y, "y", "interp1.mat", TRUE );
	write_array_matlab( xi, "xi", "interp1.mat", TRUE );
	write_array_matlab( yi, "yi", "interp1.mat", TRUE );
	array_free( x );
	array_free( y );
	array_free( xi );
	array_free( yi );
}
END_TEST

/* template
START_TEST (test_)
{
}
END_TEST
*/
Suite * init_other_suite (void){
  Suite *s = suite_create ("Other/Misc-Functions");

  TCase *tc_core = tcase_create ("OtherCore");
  tcase_add_test (tc_core, test_strip_blank );
  tcase_add_test (tc_core, test_eeg_to_stream );
  tcase_add_test (tc_core, test_nnprepare );
  tcase_add_test (tc_core, test_pqueue );
  tcase_add_test (tc_core, test_bresenham);
  tcase_add_test (tc_core, test_bresenham_segments);
  tcase_add_test (tc_core, test_disttransform);
  tcase_add_test (tc_core, test_regularization_line);
  tcase_add_test (tc_core, test_regularization_gauss);
  tcase_add_test (tc_core, test_regularization_gauss_narrowdown);
  tcase_add_test (tc_core, test_spectrogram_init);
  tcase_add_test (tc_core, test_windows);
  tcase_add_test (tc_core, test_spectgram);
  tcase_add_test (tc_core, test_check_bit);
  tcase_add_test (tc_core, test_interp1);

  tcase_set_timeout(tc_core, 20);
  suite_add_tcase (s, tc_core);
  return s;
}

