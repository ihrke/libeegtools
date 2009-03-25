/***************************************************************************
 *   Copyright (C) 2008/2009 by Matthias Ihrke   *
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

#include "nonlinear.h"


/** Get an element from phase-space reconstruction:
	 \f[
	 \vec{x}_i = \sum_{j=1}^{m}s_{i+(j-1)\tau}\vec{e}_j
	 \f]
	 \param p is the phase-space rep. of the data
	 \param i,j as in the formula above
 */
double phspace_index_ij( PhaseSpace *p, int i, int j ){
  int    idx;
  
  if( i<0 || i>p->xn || j<0 || j>=p->m ){
	 errprintf("i=%i, p->xn=%i\n, j=%i, p->m=%i\n", i, p->xn, j, p->m);
	 return 0;
  }
  idx = i - (((p->m) - (j+1)) * (p->tau));
  if( idx<0 ){
	 idx = (p->xn)+idx;
  }
  dprintf("idx=%i\n", idx);
  return(p->x[idx]);
}

/** Get an element from phase-space reconstruction:
	 \f[
	 \vec{x}_i = \sum_{j=1}^{m}s_{i+(j-1)\tau}\vec{e}_j
	 \f]
	 \param p is the phase-space rep. of the data
	 \param i as in the formula above
	 \param x output (m long vector containing the final vector)
 */
void phspace_index_i( PhaseSpace *p, int i, double *x){
  int    idx, j;

  if( i<0 || i>p->xn ){
	 errprintf("i=%i, p->xn=%i\n", i, p->xn);
	 return;
  }
  for( j=0; j<p->m; j++ ){
	 idx = i - (((p->m) - (j+1)) * (p->tau));
	 if( idx<0 ){
		idx = (p->xn)+idx;
	 }
	 x[j] = p->x[idx];
  }
}
/** Get a dimension from phase-space reconstruction:
	 \f[
	 \vec{x}_i = \sum_{j=1}^{m}s_{i+(j-1)\tau}\vec{e}_j
	 \f]
	 \param p is the phase-space rep. of the data
	 \param j as in the formula above
	 \param x output (n long vector containing the final vector)
 */
void phspace_index_j( PhaseSpace *p, int j, double *x){
  int    idx, i;

  if( j<0 || j>=p->m ){
	 errprintf("j=%i, p->m=%i\n", j, p->m);
	 return;
  }
  for( i=0; i<p->xn; i++ ){
	 idx = i - (((p->m) - (j+1)) * (p->tau));
	 if( idx<0 ){
		idx = (p->xn)+idx;
	 }
	 x[i] = p->x[idx];
  }
}

PhaseSpace* phspace_init( int m, int tau, double *x, int n ){
  PhaseSpace *p;
  p = (PhaseSpace*)malloc(sizeof(PhaseSpace));
  p->m = m;
  p->tau = tau;
  p->x = x;
  p->xn= n;
  return p;
}

void phspace_free( PhaseSpace *p ){
  free(p);
}

void phspace_print( FILE *out, PhaseSpace *p){
  fprintf( out, "PhaseSpace '%p':\n"
			  " m   = %i\n"
			  " tau = %i\n"
			  " x   = %p\n"
			  " xn  = %i\n", p, p->m, p->tau, p->x, p->xn );
}
