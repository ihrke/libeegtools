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
  dprintf("opening '%s'\n", CHECKDATADIR"/eeg061206_4_resampled500hz_filtered_DT.set");
  EEG *eeg = read_eeglab_file( CHECKDATADIR"/eeg061206_4_resampled500hz_filtered_DT.set" );
  eeg_print( stdout, eeg, 2 );

  ulonglong_t n;
  char *stream;
  stream = eeg_to_stream( eeg, NULL, &n );
  EEG *eeg2 = stream_to_eeg( stream, n, NULL );
  dprintf("short dataset: n=%llu\n", n );

  fail_unless( eeg_cmp_depth( eeg, eeg2 ) );

  eeg_print( stdout, eeg2, 2 );

  eeg_free( eeg );
  eeg_free( eeg2 );

  /* larger dataset */
  eeg = read_eeglab_file( CHECKDATADIR"/vp19_tt.set" );
  n = sizeof_eegstream( eeg );
  dprintf("long dataset: n=%llu\n", n );
  eeg_free( eeg );
}
END_TEST

START_TEST (test_nnprepare)
{
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


Suite * init_other_suite (void){
  Suite *s = suite_create ("Other/Misc-Functions");

  TCase *tc_core = tcase_create ("OtherCore");
  //  tcase_add_test (tc_core, test_strip_blank );
  //  tcase_add_test (tc_core, test_eeg_to_stream );
    tcase_add_test (tc_core, test_nnprepare );
	 //tcase_add_test (tc_core, test_pqueue );

  tcase_set_timeout(tc_core, 20);
  suite_add_tcase (s, tc_core);
  return s;
}

