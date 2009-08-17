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

#include "som.h"
#include "distances.h"


Som* som_init( int dimension, int n ){
  Som *s;
  int i;
  s=(Som*)malloc(sizeof(Som));
  s->dimension=dimension;
  s->n=n;
  s->m = (double**) malloc( n*sizeof(double*) );
  for( i=0; i<n; i++ ){
	 s->m[i] = (double*)malloc( dimension*sizeof(double));
  }
  s->d = vectordist_euclidean;
  s->nruns = 10000;
  s->initial_runs = 1000;
  s->h = neighbourhood_gaussian;
}

void   som_free( Som *s ){
  int i;
  for( i=0; i<s->n; i++ ){
	 free( s->m[i] );
  }
  free( s->m );
  free( s );
}

double neighbourhood_gaussian( int x, int m ){
  
}


