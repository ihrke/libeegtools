#include "array.h"
#include <string.h>
#include <stdio.h>


int main(){
  union u_test {
	 char c;
	 int  i;
	 float f;
	 double d;
  } utest;

  utest.d=12.4212;
  printf("%f, %i\n", utest.d, utest.f );
  /*
  int i,j,k;
  int d1=14,d2=11;
  int n=d1*d2;
  Array *a;
  float *d;
  a=array_new2( FLOAT, 2, d1, d2 );
  d = (float*)a->data;
  for( i=0; i<n; i++ ){
	 d[i] = i;
  }

  fprintf( stderr, "a=\n");
  array_print( a, -1, stderr );
  Array *b = array_slice( a, ":,:" );
  fprintf( stderr, "b=\n");
  array_print( b, -1, stderr );

  Array *c = array_slice( a, "[0 3 2 4 9 0],[1 10]");
  fprintf( stderr, "c=\n");
  array_print( c, -1, stderr );

  Array *e = array_slice( a, "[3 2 4 9 0],1-10");
  fprintf( stderr, "e=\n");
  array_print( e, -1, stderr );

  int x;
  Dtype dt=INT;
  x = (array_DTYPE_TO_CDTYPE( INT ))2.3;
  printf("x=%i\n", x );

  array_free( a );
  array_free( b );
  
  */


  return 0;
}
