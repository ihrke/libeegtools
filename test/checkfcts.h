/*
** checkfcts.h
** 
** Made by Matthias Ihrke
** Login   <mihrke@localhost>
** 
** Started on  Thu Oct 11 13:18:06 2007 Matthias Ihrke
** Last update Sun Oct 14 14:27:05 2007 Matthias Ihrke
*/

#ifndef   	CHECKFCTS_H_
# define   	CHECKFCTS_H_

#define CHECKDATADIR "./data"
#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "helper.h"

int cmpdouble(double d1, double d2, int precision);
int isequal_doublearray_binfile(double *d, int n, const char *filename, int precision);
int isequal_doublearray(const double *d1, const double *d2, int n, int precision);

#endif 	    /* !CHECKFCTS_H_ */
