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

int cmpdouble(double d1, double d2, int precision){
  double epsilon;
  epsilon = pow(10, -1*precision);
 /*  printf("epsilon=%.10f, d1=%.10f, d2=%.10f\n", epsilon, d1, d2); */
/*   printf("fabs(d1-d2)=%.10f, e*fabs(d1)=%.10f\n", fabs(d1 - d2), epsilon * fabs(d1)); */
  if(fabs(d1 - d2) < epsilon * fabs(d1))
    return 0;
  else if(d1<d2) return -1;
  else return 1;
}

int isequal_doublearray(const double *d1, const double *d2, int n, int precision){

}

int isequal_doublearray_binfile(double *d, int n, const char *filename, int precision){
  FILE *f;
  double *d2;
  int i, flag=1;

  d2 = (double*)malloc(n*sizeof(double));
  if((f=fopen(filename, "rb"))==NULL) return -1;
  if(fread(d, n*sizeof(double), 1, f)<n) return -1;
  fclose(f);

  for(i=0; i<n; i++){
    printf("d[%i]=%.10f, d2[%i]=%.10f\n", i, d[i], i, d2[i]);
    if(cmpdouble(d[i], d2[i], precision)!=0){
      flag=0; break;
    }
  }
  free(d2);
  return flag;
}

