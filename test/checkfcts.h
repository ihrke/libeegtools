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
#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/** Compare two doubles to a certain precision.
 * cmpdouble()
 *
 * \param precision to which position after comma it is compared
 * \return 
 *   - -1 if d1<d2
 *   - 0 if d1==d2
 *   - 1 if d1>d2
 */
int cmpdouble(double d1, double d2, int precision);

/** Compare an array of doubles to doubles stored binary file.
 * The data must be stored in 32 bit floats
 *
 * \param precision to which position after comma it is compared
 * \return 
 *   - 0 if there is an element that is not equal
 *   - 1 if everything is equal
 * \todo check it, does not seem to work properly
 */
int isequal_doublearray_binfile(double *d, int n, const char *filename, int precision);

/** Compares two double arrays.
 * isequal_doublearray()
 *
 * \param precision \see cmpdouble()
 * \return
 *   - 0 equal they are
 *   - 1 otherwhise
 * \todo implement it
 */
int isequal_doublearray(const double *d1, const double *d2, int n, int precision);

#endif 	    /* !CHECKFCTS_H_ */
