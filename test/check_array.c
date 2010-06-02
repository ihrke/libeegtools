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
#include "array.h"
#include "linalg.h"
#include "helper.h"

START_TEST (test_sizeof)
{
  int s;
  array_SIZEOF_DTYPE( s, CHAR );
  fail_unless( s==sizeof(char) );
  array_SIZEOF_DTYPE( s, INT );
  fail_unless( s==sizeof(int) );
  array_SIZEOF_DTYPE( s, LONG );
  fail_unless( s==sizeof(long) );
  array_SIZEOF_DTYPE( s, ULONG );
  fail_unless( s==sizeof(unsigned long) );
  array_SIZEOF_DTYPE( s, FLOAT );
  fail_unless( s==sizeof(float) );
  array_SIZEOF_DTYPE( s, DOUBLE );
  fail_unless( s==sizeof(double) );
}
END_TEST

START_TEST (test_dtstring)
{
  char *s;
  array_DTYPESTRING( s, CHAR );
  fail_unless( strcmp( s, "char" )==0 );
  array_DTYPESTRING( s, INT );
  fail_unless( strcmp( s, "int" )==0 );
  array_DTYPESTRING( s, LONG );
  fail_unless( strcmp( s, "long" )==0 );
  array_DTYPESTRING( s, ULONG );
  fail_unless( strcmp( s, "unsigned long" )==0 );
  array_DTYPESTRING( s, FLOAT );
  fail_unless( strcmp( s, "float" )==0 );
  array_DTYPESTRING( s, DOUBLE );
  fail_unless( strcmp( s, "double" )==0 );
}
END_TEST



     
START_TEST (test_array_new)
{
  int size[2]={10,11};
  Array *a=array_new( DOUBLE, 2, size );
  fail_unless( a->dtype==DOUBLE );
  fail_unless( a->dtype_size==sizeof(double) );
  fail_unless( a->ndim==2 );
  fail_unless( a->size[0]==10 );
  fail_unless( a->size[1]==11 );
  fail_unless( a->nbytes==10*11*sizeof(double) );
  array_free( a );
}
END_TEST

START_TEST (test_array_new2)
{
  Array *a=array_new2( DOUBLE, 2, 10, 11 );
  fail_unless( a->dtype==DOUBLE );
  fail_unless( a->dtype_size==sizeof(double) );
  fail_unless( a->ndim==2 );
  fail_unless( a->size[0]==10 );
  fail_unless( a->size[1]==11 );
  fail_unless( a->nbytes==10*11*sizeof(double) );
  array_free( a );
}
END_TEST

START_TEST (test_array_fromptr)
{
  double *m;
  int i,n=100;
  Array *a;
  MALLOC( m, n, double );
  for( i=0; i<n; i++ ){
	 m[i] = i;
  }
  a=array_fromptr2( DOUBLE, 2, m, 10, 10 );

  fail_unless( a->dtype==DOUBLE );
  fail_unless( a->dtype_size==sizeof(double) );
  fail_unless( a->ndim==2 );
  fail_unless( a->size[0]==10 );
  fail_unless( a->size[1]==10 );
  fail_unless( a->nbytes==10*10*sizeof(double) );
  array_free( a );
}
END_TEST

START_TEST (test_index1)
{
  int i,n=110;
  Array *a;
  double *d;
  a=array_new2( DOUBLE, 1, n );
  d = (double*)a->data;
  for( i=0; i<n; i++ ){
	 d[i] = i;
  }
  for( i=0; i<n; i++ ){
	 fail_unless( d[i]==array_INDEX1(a, double, i) );
  }
  array_free( a );
}
END_TEST

START_TEST (test_index2)
{
  int i,j,k;
  int d1=14,d2=11;
  int n=d1*d2;
  Array *a;
  int *d;
  a=array_new2( INT, 2, d1, d2 );
  d = (int*)a->data;
  for( i=0; i<n; i++ ){
	 d[i] = i;
  }
  k=0;
  for( i=0; i<d1; i++ ){
	 for( j=0; j<d2; j++ ){
		fail_unless( d[k]==array_INDEX2(a, int, i, j),
						 "%i!=%i",d[k], array_INDEX2(a, int, i, j) );
		k++;
	 }
  }
  array_free( a );
}
END_TEST


START_TEST (test_index3)
{
  int i,j,k,l; 
  int d1=10, d2=11, d3=12;
  int n=d1*d2*d3;

  Array *a;
  int *d;
  a=array_new2( INT, 3, d1, d2, d3 );
  d = (int*)a->data;
  for( i=0; i<n; i++ ){
	 d[i] = i;
  }
  l=0;
  for( i=0; i<d1; i++ ){
	 for( j=0; j<d2; j++ ){
		for( k=0; k<d3; k++ ){
		  fail_unless( d[l]==array_INDEX3(a, int, i, j, k),
							"%i!=%i, idx=%i\n", d[l],array_INDEX3(a, int, i, j, k), l );
		  l++;
		}
	 }
  }
  array_free( a );
}
END_TEST

START_TEST (test_index)
{
  int i,j,k,l,m; 
  int d1=10, d2=11, d3=12, d4=13;
  int n=d1*d2*d3*d4;

  Array *a;
  int *d;
  a=array_new2( INT, 4, d1, d2, d3, d4 );
  d = (int*)a->data;
  for( i=0; i<n; i++ ){
	 d[i] = i;
  }
  m=0;
  for( i=0; i<d1; i++ ){
	 for( j=0; j<d2; j++ ){
		for( k=0; k<d3; k++ ){
		  for( l=0; l<d4; l++ ){
			 fail_unless( d[m]==*(int*)array_index2(a, i, j, k, l),
							  "%i!=%i, idx=%i\n", 
							  d[m],*(int*)array_index2(a,i, j, k, l), m );
			 m++;
		  }
		}
	 }
  }
  array_free( a );
}
END_TEST

START_TEST (test_index_dbl)
{
  ulong i,j,k,l; 
  uint d1=10, d2=11, d3=12, d4=13;
  ulong n=d1*d2*d3*d4;
  ulong m;
  double r;
  Array *a;
  double *d;
  a=array_new2( DOUBLE, 4, d1, d2, d3, d4 );
  d = (double*)a->data;
  for( i=0; i<n; i++ ){
	 d[i] = i;
  }
  m=0;
  for( i=0; i<d1; i++ ){
	 for( j=0; j<d2; j++ ){
		for( k=0; k<d3; k++ ){
		  for( l=0; l<d4; l++ ){
			 r = *(double*)array_index2(a, i, j, k, l);
			 fail_unless( cmpdouble(d[m],r,2)==0,
							  "%f!=%f, idx=%li/%li\n", 
							  d[m],r, m,n );
			 m++;
		  }
		}
	 }
  }
  array_free( a );
}
END_TEST

START_TEST (test_slice)
{
  int i,j,k;
  int d1=14,d2=11;
  int n=d1*d2;
  Array *a;
  int *d;
  a=array_new2( INT, 2, d1, d2 );
  d = (int*)a->data;
  for( i=0; i<n; i++ ){
	 d[i] = i;
  }
  
  Array *b = array_slice( a, ":,:" );
  fail_unless( memcmp( a->data, b->data, a->nbytes )==0 );

  Array *c = array_slice( a, "[0 3 2 4 9 0],[1 10]");
  fail_unless( c->ndim==2 );
  fail_unless( c->size[0]==6 && c->size[1]==2 );
  Array *e = array_slice( a, "[3 2 4 9 0],1-10");
  fail_unless( e->ndim==2 );
  fail_unless( e->size[0]==5 && e->size[1]==10 );

  Array *f = array_slice( a, "3,1-10");
  fail_unless( f->ndim==1 );
  fail_unless( f->size[0]==10 );

  array_free( a );
  array_free( b );
  array_free( c );
  array_free( e );
  array_free( f );
}
END_TEST

START_TEST (test_scale)
{
  Array *a=array_new_dummy( INT, 3, 3,4,5 );
  long i;
  array_scale( a, 2 );
  for( i=0; i<a->nbytes/a->dtype_size; i++ ){
	 fail_unless( *((int*)(a->data+i*a->dtype_size))==(int)i*2,
					  "%i!=%i", *((int*)a->data+i*a->dtype_size),(int)i*2 );
  }
  Array *b=array_new_dummy( DOUBLE, 3, 3,4,5 );
  array_scale( b, 0.5 );
  for( i=0; i<b->nbytes/b->dtype_size; i++ ){
	 fail_unless( cmpdouble( *((double*)(b->data+i*b->dtype_size)), (double)i*0.5, 3 )==0,
					  "%.2f!=%.2f",*((double*)(b->data+i*b->dtype_size)), (double)i*0.5 );
  }
  array_free( a );
  array_free( b );
}
END_TEST

START_TEST (test_copy)
{
  Array *a=array_new_dummy( FLOAT, 3, 3,4,5 );
  Array *b=array_copy( a, TRUE );
  Array *c=array_copy( a, FALSE );

  fail_unless( b->data != a->data );
  fail_unless( cmpdouble( array_INDEX3( b, float, 0,0,0),
								  array_INDEX3( a, float, 0,0,0), 3)==0);
  fail_unless( b->free_data==TRUE );
  fail_unless( c->data == a->data );
  fail_unless( cmpdouble( array_INDEX3( c, float, 0,0,0),
								  array_INDEX3( a, float, 0,0,0), 3)==0);
  fail_unless( c->free_data==FALSE );

  array_free( a );
  array_free( b );
  array_free( c );
}
END_TEST

START_TEST (test_comparable)
{
  Array *a=array_new_dummy( FLOAT, 3, 3,4,5 );
  Array *b=array_copy( a, TRUE );
  fail_unless( array_comparable( a, b ) );
  b->size[0]=2;
  fail_if( array_comparable( a, b ) );
  array_free( a );
  array_free( b );
}
END_TEST

START_TEST (test_reverse)
{
  Array *a=array_new_dummy( DOUBLE, 1, 10 );
  Array *b=array_new_dummy( DOUBLE, 1, 11 );
  Array *c=array_copy( a, TRUE );
  Array *d=array_copy( b, TRUE );
  Array *e=array_new_dummy( DOUBLE, 3, 10, 11, 12 );
  Array *f=array_copy( e, TRUE );
  
  array_reverse( a );
  array_reverse( b );
  array_reverse( f );

  int i;
  for( i=0; i<10; i++ ){
	 fail_unless( cmpdouble( array_INDEX1( c, double, i ),
									 array_INDEX1( a, double, 10-i-1), 2 )==0 );
  }
  for( i=0; i<11; i++ ){
	 fail_if( cmpdouble( array_INDEX1( d, double, i ),
								array_INDEX1( b, double, 11-i-1), 3)!=0 );
  }

  for( i=0; i<array_NUMEL( e ); i++ ){
	 fail_unless( cmpdouble( array_INDEX1( e, double, i ),
									 array_INDEX1( f, double, array_NUMEL(e)-i-1), 2 )==0 );	 
  }

  array_free( a );
  array_free( b );
  array_free( c );
  array_free( d );
  array_free( e );
  array_free( f );
}
END_TEST

START_TEST (test_random)
{
  Array *a=array_randunif( 0, 1, 100 );
  int i;
  for( i=0; i<array_NUMEL( a ); i++ ){
	 fail_if( array_INDEX1( a, double, i )<0 ||
				 array_INDEX1( a, double, i )>1 );
  }
  Array *b=array_randunif( 0, 3, 10, 11, 12 );
  for( i=0; i<array_NUMEL( b ); i++ ){
	 fail_if( array_INDEX1( b, double, i )<0 ||
				 array_INDEX1( b, double, i )>1 );
  }

  array_free(a);
  array_free(b);
}
END_TEST

START_TEST (test_shuffle)
{
  int i;
  Array *a=array_new_dummy( DOUBLE, 1, 100 );
  Array *b=array_new_dummy( UINT, 3, 11, 12, 13 );
  array_shuffle( a, 0 );
  array_shuffle( b, 0 );

  for( i=0; i<array_NUMEL(a); i++ ){
	 fail_if( array_INDEX1( a, double, i )<0 ||
				 array_INDEX1( a, double, i )>=array_NUMEL(a) );
  }
  for( i=0; i<array_NUMEL(b); i++ ){
	 fail_if( array_INDEX1( b, double, i )<0 ||
				 array_INDEX1( b, double, i )>=array_NUMEL(b) );
  }

  array_free( a );
  array_free( b );
}
END_TEST

START_TEST (test_rowindex)
{
  int d1=2, d2=3, d3=4;
  Array *a=array_new_dummy( UINT, 3, d1, d2, d3 );

  int i,j,k,offset=0;
  int index[3];
  for( i=0; i<d1; i++ ){
	 for( j=0; j<d2; j++ ){
		for( k=0; k<d3; k++ ){
		  array_calc_rowindex( offset, a->size, a->ndim, index );
		  /* dprintf("i,j,k=(%i,%i,%i), index=(%i,%i,%i), offset=%i\n", */
		  /* 			 i,j,k,index[0],index[1],index[2],offset); */
		  fail_unless( index[0]==i );
		  fail_unless( index[1]==j );
		  fail_unless( index[2]==k );
		  offset++;
		}
	 }
  }

  array_free( a );
}
END_TEST

START_TEST (test_colindex)
{
  int d1=2, d2=3, d3=4;
  Array *a=array_new_dummy( UINT, 3, d1, d2, d3 );

  int i,j,k,offset=0;
  int index[3];
  for( i=0; i<d3; i++ ){
	 for( j=0; j<d2; j++ ){
		for( k=0; k<d1; k++ ){
		  array_calc_colindex( offset, a->size, a->ndim, index );
		  /* dprintf("i,j,k=(%i,%i,%i), index=(%i,%i,%i), offset=%i\n", */
		  /* 			 i,j,k,index[0],index[1],index[2],offset); */
		  fail_unless( index[0]==k );
		  fail_unless( index[1]==j );
		  fail_unless( index[2]==i );
		  offset++;
		}
	 }
  }

  array_free( a );
}
END_TEST
START_TEST (test_convertcolrow)
{
  int d1=2, d2=3, d3=4;
  Array *a=array_new_dummy( UINT, 3, d1, d2, d3 );
  Array *b=array_convert_rowcolmajor( a, TRUE );

  int i,j,k;
  int offset=0;
  int index[3];
  for( i=0; i<d3; i++ ){
	 for( j=0; j<d2; j++ ){
		for( k=0; k<d1; k++ ){
		  dprintf("a(%i,%i,%i)=%i, a+(offset=%i)=%i, b+offset=%i\n",
					 i,j,k,  array_INDEX3( a, uint,i,j,k),
					 offset, array_INDEX1( a,uint,offset),
					 array_INDEX1( b,uint,offset) );
		  fail_unless( array_INDEX1( b,uint,offset) ==
							array_INDEX3( a, uint,k,j,i ) );
			 
		  /* array_calc_colindex( array_INDEX1( b,uint,offset), */
		  /* 							  a->size, a->ndim, index ); */
		  
		  /* fail_unless( index[0]==i, "index[0]=%i, i=%i", index[0], i ); */
		  /* fail_unless( index[1]==j, "index[1]=%i, j=%i", index[1], j ); */
		  /* fail_unless( index[2]==k, "index[2]=%i, k=%i", index[2], k ); */

		  offset++;
		}
	 }
  }

  array_free( a );
  array_free( b );
}
END_TEST

START_TEST (test_concatenate_rows)
{
  Array *a=array_new_dummy( UINT, 2, 10, 3 );
  Array *b=array_new_dummy( UINT, 2, 4, 3 );
  Array *c=array_concatenate( a, b, 0 );

  
  /* array_print( a, -1, stderr ); */
  /* array_print( b, -1, stderr ); */
  /* array_print( c, -1, stderr ); */

  fail_unless( c!=NULL );
  int i,j;
  for( i=0; i<a->size[0]; i++ ){
	 for( j=0; j<a->size[1]; j++ ){
		fail_unless( array_INDEX2( c, uint,i,j )==array_INDEX2( a, uint,i,j ) );
	 }
  }
  for( i=0; i<b->size[0]; i++ ){
	 for( j=0; j<b->size[1]; j++ ){
		fail_unless( array_INDEX2( c, uint,i+a->size[0],j )==array_INDEX2( b, uint,i,j ) );
	 }
  }

  array_free( a );
  array_free( b );
  array_free( c );
}
END_TEST

START_TEST (test_concatenate_cols)
{
  Array *a=array_new_dummy( UINT, 2, 10, 3 );
  Array *b=array_new_dummy( UINT, 2, 10, 4 );
  Array *c=array_concatenate( a, b, 1 );
  
  /* array_print( a, -1, stderr ); */
  /* array_print( b, -1, stderr ); */
  /* array_print( c, -1, stderr ); */

  fail_unless( c!=NULL );
  int i,j;
  for( i=0; i<a->size[0]; i++ ){
	 for( j=0; j<a->size[1]; j++ ){
		fail_unless( array_INDEX2( c, uint,i,j )==array_INDEX2( a, uint,i,j ) );
	 }
  }
  for( i=0; i<b->size[0]; i++ ){
	 for( j=0; j<b->size[1]; j++ ){
		fail_unless( array_INDEX2( c, uint,i,j+a->size[1] )==array_INDEX2( b, uint,i,j ) );
	 }
  }

  array_free( a );
  array_free( b );
  array_free( c );
}
END_TEST

START_TEST (test_concatenate_cols_dbl)
{
  Array *a=array_new_dummy( DOUBLE, 2, 10, 3 );
  Array *b=array_new_dummy( DOUBLE, 2, 10, 4 );
  Array *c=array_concatenate( a, b, 1 );
  
  /* array_print( a, -1, stderr );  */
  /* array_print( b, -1, stderr ); */
  /* array_print( c, -1, stderr );  */

  fail_unless( c!=NULL );
  int i,j;
  for( i=0; i<a->size[0]; i++ ){
	 for( j=0; j<a->size[1]; j++ ){
		fail_unless( cmpdouble( array_INDEX2( c, uint,i,j ),
										array_INDEX2( a, uint,i,j ),3)==0 );
	 }
  }
  for( i=0; i<b->size[0]; i++ ){
	 for( j=0; j<b->size[1]; j++ ){
		fail_unless( cmpdouble( array_INDEX2( c, uint,i,j+a->size[1] ),
										array_INDEX2( b, uint,i,j ),3)==0 );
	 }
  }

  Array *d=array_concatenate( a, NULL, DIM_COLS );
  fail_unless( array_comparable( a, d ) );
  for( i=0; i<array_NUMEL(a); i++ ){
	 fail_unless( cmpdouble( array_INDEX1( a,double,i),
									 array_INDEX1( d,double,i), 4)==0 );
  }

  array_free( a );
  array_free( b );
  array_free( c );
  array_free( d );
}
END_TEST

START_TEST (test_min)
{
  Array *a=array_new_dummy( DOUBLE, 2, 10, 3 );
  double m=*(double*)array_min( a );
  fail_unless( cmpdouble( m, 0, 4 )==0 );

  Array *b=array_new_dummy( UINT, 2, 10, 3 );
  uint m2=*(uint*)array_min( b );
  fail_unless( m2==0 );

  mat_IDX( a, 2, 4 )=-10.0;
  m=*(double*)array_min( a );
  fail_unless( cmpdouble( m, -10, 4)==0 );

  array_free( a );
  array_free( b );
}
END_TEST

START_TEST (test_max)
{
  Array *a=array_new_dummy( DOUBLE, 2, 10, 3 );
  double m=*(double*)array_max( a );
  fail_unless( cmpdouble( m, (double)10*3-1, 4 )==0 );

  Array *b=array_new_dummy( UINT, 2, 10, 3 );
  uint m2=*(uint*)array_max( b );
  fail_unless( m2==10*3-1 );

  mat_IDX( a, 2, 4 )=100.0;
  m=*(double*)array_max( a );
  fail_unless( cmpdouble( m, 100, 4)==0 );

  array_free( a );
  array_free( b );
}
END_TEST


/* template
START_TEST (test_)
{
}
END_TEST
*/

Suite * init_array_suite (void){
  Suite *s = suite_create ("Numerical Array-Functions");

  TCase *tc_core = tcase_create ("ArrayCore");
  tcase_add_test (tc_core, test_sizeof );
  tcase_add_test (tc_core, test_dtstring );
  tcase_add_test (tc_core, test_array_new );
  tcase_add_test (tc_core, test_array_new2 );
  tcase_add_test (tc_core, test_array_fromptr );
  tcase_add_test (tc_core, test_index  );
  tcase_add_test (tc_core, test_index_dbl  );
  tcase_add_test (tc_core, test_index1 );
  tcase_add_test (tc_core, test_index2 );
  tcase_add_test (tc_core, test_index3 );
  tcase_add_test (tc_core, test_slice );
  tcase_add_test (tc_core, test_scale );
  tcase_add_test (tc_core, test_copy);
  tcase_add_test (tc_core, test_comparable );
  tcase_add_test (tc_core, test_reverse );
  tcase_add_test (tc_core, test_random);
  tcase_add_test (tc_core, test_shuffle);
  tcase_add_test (tc_core, test_rowindex);
  tcase_add_test (tc_core, test_colindex);
  tcase_add_test (tc_core, test_convertcolrow);
  tcase_add_test (tc_core, test_concatenate_rows);
  tcase_add_test (tc_core, test_concatenate_cols);
  tcase_add_test (tc_core, test_concatenate_cols_dbl);
  tcase_add_test (tc_core, test_min);
  tcase_add_test (tc_core, test_max);

  tcase_set_timeout(tc_core, 20);
  suite_add_tcase (s, tc_core);
  return s;
}

