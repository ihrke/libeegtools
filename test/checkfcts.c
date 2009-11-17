/**
* \file checkfcts.c
* \brief functions that help to check other function.
**
\code 
** Made by (Matthias Ihrke)
** Login   <mihrke@localhost>
** 
** Started on  Thu Oct 11 13:17:59 2007 Matthias Ihrke
** Last update Sun Oct 14 14:25:21 2007 Matthias Ihrke
\endcode
*/

#include "checkfcts.h"

/** Compare two doubles to a certain precision.
 * cmpdouble()
 *
 * \param precision to which position after comma it is compared
 * \return 
 *   - -1 if d1<d2
 *   - 0 if d1==d2
 *   - 1 if d1>d2
 */
int cmpdouble(double d1, double d2, int precision){
  double epsilon;
  epsilon = pow(10, -1*precision);
 /*  printf("epsilon=%.10f, d1=%.10f, d2=%.10f\n", epsilon, d1, d2); */
/*   printf("fabs(d1-d2)=%.10f, e*fabs(d1)=%.10f\n", fabs(d1 - d2), epsilon * fabs(d1)); */
//  dprintf("%f==%f with epsilon=%f, d1-d2=%f, |d1-d2|=%f, %i?\n", d1, d2, epsilon, d1-d2, fabs(d1-d2), fabs(d1-d2)<epsilon );
  if(fabs(d1 - d2) < epsilon*fabs(d1) ){
    return 0;
  } else if(d1<d2){
	 return -1;
  } else {
	 return 1;
  }
}
/** Compares two double arrays.
 * isequal_doublearray()
 *
 * \param precision \see cmpdouble()
 * \return
 *   - 0 equal they are
 *   - 1 otherwhise
 * \todo implement it
 */
int isequal_doublearray(const double *d1, const double *d2, int n, int precision){
  int i, res;
  for( i=0; i<n; i++ ){
	 res=cmpdouble(d1[i], d2[i], precision);
	 if( res )
		return res;
  }
  return 0;
}

/** Compare an array of doubles to doubles stored binary file.
 * The data must be stored in 32 bit floats
 *
 * \param precision to which position after comma it is compared
 * \return 
 *   - 0 if everything is equal
 *   - 1 if there is an element that is not equal
 * \todo check it, does not seem to work properly
 */
int isequal_doublearray_binfile(double *d, int n, const char *filename, int precision){
  FILE *f;
  double *d2;
  int i, flag=1;

  d2 = (double*)malloc(n*sizeof(double));
  if((f=fopen(filename, "rb"))==NULL){
	 errprintf("opening file '%s' failed\n", filename);
	 return -1;
  }

  ffread(d2, sizeof(double), n, f);
  fclose(f);

  flag = isequal_doublearray( d, d2, n, precision );
  free(d2);
  return flag;
}

