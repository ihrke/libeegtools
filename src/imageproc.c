/***************************************************************************
 *   Copyright (C) 2010 by Matthias Ihrke   *
 *   ihrke@nld.ds.mpg.de
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
#include <float.h>
#include "helper.h"
#include "mathadd.h"
#include "imageproc.h"
#include "linalg.h"

/** \brief calculate the (signed) distance-transform of a binary image.

	 This distance-transform (DT) of a binary image includes at each
	 pixel the shortest distance to the closest nonzero point in
	 the input image.

	 Example (left is the binary matrix, right the DT; in the upper right corner is 
	 another white pixel):
	 \image html disttransform.jpg

	 This implementation uses the dead-reckoning algorithm from
	 \verbatim
	 @article{grevera2004dead,
  	   title={{The" Dead reckoning" signed distance transform}},
	   author={Grevera, G.J.},
	   journal={Computer Vision and Image Understanding},
	   volume={95},
	   number={3},
	   pages={317--333},
	   year={2004}
	 }
	 \endverbatim

	 \param in INT-array representing a binary image
	 \param dt DOUBLE-array of same dimensions as in, or NULL (alloc in function)
	 \return ptr to the distance transform of in
 */
Array*  disttransform_deadreckoning( const Array *in, Array *dt ){ 
  int x, y;
  int X, Y;
  double d1=1;
  double d2=sqrt(2);
 
  if( in->dtype!=INT || in->ndim!=2 ){
	 errprintf("Need a 2-dim INT array as input\n");
	 return NULL;
  }
  X = in->size[0];
  Y = in->size[1];
  dprintf("X,Y=%i,%i\n",X,Y);
  if( !dt ){
	 dt=array_new( DOUBLE, 2, in->size );
  }

  /* tmp arrays */
  Array *P1 = array_new( DOUBLE, 2, in->size );
  Array *P2 = array_new( DOUBLE, 2, in->size );

  /* init */
  for(y=0; y<Y; y++){
	 for(x=0; x<X; x++){
		mat_IDX(dt,x,y) = DBL_MAX;
		mat_IDX(P1,x,y)=-1;
		mat_IDX(P2,x,y)=-1;
		if(x>0 && y>0 && y<Y-1 && x<X-1){
		  if( array_INDEX2(in,int,x-1,y)!=array_INDEX2(in,int,x,y) || 
				array_INDEX2(in,int,x+1,y)!=array_INDEX2(in,int,x,y) ||
				array_INDEX2(in,int,x,y-1)!=array_INDEX2(in,int,x,y) || 
				array_INDEX2(in,int,x,y+1)!=array_INDEX2(in,int,x,y) ){
			 mat_IDX(dt,x,y)=0;
			 mat_IDX(P1,x,y)=x;
			 mat_IDX(P2,x,y)=y;
		  }
		}            
	 }
  }
  /*-------------------- computation -----------------------*/
  /* first pass */
  for(y=1; y<Y-1; y++){
	 for(x=1; x<X-1; x++){
		if( mat_IDX(dt,x-1,y-1)+d2 < mat_IDX(dt,x,y) ){
		  mat_IDX(P1,x,y) = mat_IDX(P1,x-1,y-1);
		  mat_IDX(P2,x,y) = mat_IDX(P2,x-1,y-1);
		  mat_IDX(dt,x,y) = sqrt( SQR( x-mat_IDX(P1,x,y) ) + 
										  SQR( y-mat_IDX(P2,x,y) ) );
		}
		if( mat_IDX(dt,x,y-1)+d1 < mat_IDX(dt,x,y) ) {
		  mat_IDX(P1,x,y) = mat_IDX(P1,x,y-1);
		  mat_IDX(P2,x,y) = mat_IDX(P2,x,y-1);
		  mat_IDX(dt,x,y) = sqrt( SQR( x-mat_IDX(P1,x,y) ) + 
										  SQR( y-mat_IDX(P2,x,y) ) );		 
		}    
		if( mat_IDX(dt,x+1,y-1)+d2 < mat_IDX(dt,x,y) ) {
		  mat_IDX(P1,x,y) = mat_IDX(P1,x+1,y-1);
		  mat_IDX(P2,x,y) = mat_IDX(P2,x+1,y-1);
		  mat_IDX(dt,x,y) = sqrt( SQR( x-mat_IDX(P1,x,y) ) + 
										  SQR( y-mat_IDX(P2,x,y) ) );		 
		}    
		if( mat_IDX(dt,x-1,y)+d1 < mat_IDX(dt,x,y) ) {
		  mat_IDX(P1,x,y) = mat_IDX(P1,x-1,y);
		  mat_IDX(P2,x,y) = mat_IDX(P2,x-1,y);
		  mat_IDX(dt,x,y) = sqrt( SQR( x-mat_IDX(P1,x,y) ) + 
										  SQR( y-mat_IDX(P2,x,y) ) );		 
		}    
	 }
  }
  
  
  /* final pass */
  for(y=Y-2; y>0; y--){
	 for(x=X-2; x>0; x--){
		if( mat_IDX(dt,x+1,y)+d1 < mat_IDX(dt,x,y) ){
		  mat_IDX(P1,x,y) = mat_IDX(P1,x+1,y);
		  mat_IDX(P2,x,y) = mat_IDX(P2,x+1,y);
		  mat_IDX(dt,x,y) = sqrt( SQR( x-mat_IDX(P1,x,y) ) + 
										  SQR( y-mat_IDX(P2,x,y) ) );
		}
		if( mat_IDX(dt,x-1,y+1)+d2 < mat_IDX(dt,x,y) ){
		  mat_IDX(P1,x,y) = mat_IDX(P1,x-1,y+1);
		  mat_IDX(P2,x,y) = mat_IDX(P2,x-1,y+1);
		  mat_IDX(dt,x,y) = sqrt( SQR( x-mat_IDX(P1,x,y) ) + 
										  SQR( y-mat_IDX(P2,x,y) ) );
		}
		if( mat_IDX(dt,x,y+1)+d1 < mat_IDX(dt,x,y) ){
		  mat_IDX(P1,x,y) = mat_IDX(P1,x,y+1);
		  mat_IDX(P2,x,y) = mat_IDX(P2,x,y+1);
		  mat_IDX(dt,x,y) = sqrt( SQR( x-mat_IDX(P1,x,y) ) + 
										  SQR( y-mat_IDX(P2,x,y) ) );
		}
		if( mat_IDX(dt,x+1,y+1)+d2 < mat_IDX(dt,x,y) ){
		  mat_IDX(P1,x,y) = mat_IDX(P1,x+1,y+1);
		  mat_IDX(P2,x,y) = mat_IDX(P2,x+1,y+1);
		  mat_IDX(dt,x,y) = sqrt( SQR( x-mat_IDX(P1,x,y) ) + 
										  SQR( y-mat_IDX(P2,x,y) ) );
		}
	 }
  }
  

  /* indicate in- and outside */
  for(y=Y-2; y>0; y--){
	 for(x=X-2; x>0; x--){
		if( array_INDEX2(in,int,x,y)>0 ) {
		  mat_IDX(dt,x,y) = -1*mat_IDX(dt,x,y);
		}
	 }
  }  

  /* corners */
  mat_IDX(dt,0,  0  )=mat_IDX(dt,1,  1);
  mat_IDX(dt,0,  Y-1)=mat_IDX(dt,1,  Y-2);
  mat_IDX(dt,X-1,Y-1)=mat_IDX(dt,X-2,Y-2);
  mat_IDX(dt,X-1,0  )=mat_IDX(dt,X-2,1);

  /* borders */
  for(x=1; x<X-1; x++){
  	 mat_IDX(dt,x,0)   = mat_IDX(dt,x,1);
  	 mat_IDX(dt,x,Y-1) = mat_IDX(dt,x,Y-2);
  }
  for(y=1; y<Y-1; y++){
  	 mat_IDX(dt,0,y)   = mat_IDX(dt,1,y);
  	 mat_IDX(dt,X-1,y) = mat_IDX(dt,X-2,y);
  }
  /*-------------------- /computation -----------------------*/

  array_free( P1 );
  array_free( P2 );

  return dt;
}

/** \brief Connected line-segments with bresenhams Line-drawing-algorithm.

	 Calls bresenham_line() for each pair of points.

	 Example:

	 \image html bresenham_segments.jpg
	 
	 \param  points 2xM INT-array containig coordinates points through which the line-segments
	      should lead
	 \return a pointer to the line saved in pairs (x,y) of the coordinates (2 x N INT-Array);
*/
Array*  bresenham_linesegments( const Array *points ){
  if( points->dtype!=INT ){
	 errprintf("Need INT-array for coordinates\n"); return NULL;
  }
  if( points->ndim!=2 ){
	 errprintf("Need 2D-array for coordinates\n"); return NULL;
  }
  if( points->size[1]<2 ){
	 errprintf("Need 2 coordinates per point\n"); return NULL;
  }
  int start[2], end[2];
  Array *bres=NULL, *line, *tmp;
  int i;
  for( i=0; i<points->size[1]-1; i++ ){
	 start[0]=array_INDEX2( points,int,0,i );
	 start[1]=array_INDEX2( points,int,1,i );
	 end[0]=array_INDEX2( points,int,0,i+1 );
	 end[1]=array_INDEX2( points,int,1,i+1 );

	 dprintf("(%i,%i)->(%i,%i)\n",
				start[0],start[1],end[0],end[1]);
	 line=bresenham_line( start, end );
	 tmp=bres;
	 bres=array_concatenate( tmp, line, DIM_COLS );
	 array_free( tmp );
	 array_free( line );
  }

  return bres;
}

/** \brief Bresenhams Line-drawing-algorithm.

	 (stolen and modified from Wikipedia)
	 
	 \param  start  coordinates of starting point
	 \param  end    coordinates of last point
	 \return a pointer to the line saved in pairs (x,y) of the coordinates (2 x N INT-Array);
*/
Array*  bresenham_line( int start[2], int end[2] ){
  int x, y, t, dx, dy, incx, incy, pdx, pdy, ddx, ddy, es, el, err;

  /* Entfernung in beiden Dimensionen berechnen */
  dx = end[0] - start[0];
  dy = end[1] - start[1];
  
  /* Vorzeichen des Inkrements bestimmen */
  incx = sgn(dx);
  incy = sgn(dy);
  if(dx<0) dx = -dx;
  if(dy<0) dy = -dy;

   /* feststellen, welche Entfernung größer ist */
  if (dx>dy)
	 {
      /* x ist schnelle Richtung */
      pdx=incx; pdy=0;    /* pd. ist Parallelschritt */
      ddx=incx; ddy=incy; /* dd. ist Diagonalschritt */
      es =dy;   el =dx;   /* Fehlerschritte schnell, langsam */
	 } else
	 {
      /* y ist schnelle Richtung */
      pdx=0;    pdy=incy; /* pd. ist Parallelschritt */
      ddx=incx; ddy=incy; /* dd. ist Diagonalschritt */
      es =dx;   el =dy;   /* Fehlerschritte schnell, langsam */
	 }
 
  /* Initialisierungen vor Schleifenbeginn */
  x = start[0];
  y = start[1];
  err = el/2;

  Array *out = array_new2( INT, 2, 2, el+1 );

  array_INDEX2(out,int,0,0)=x;
  array_INDEX2(out,int,1,0)=y;
  /*------------------ computation -------------------------*/
  for(t=1; t<=el; t++) {
	 /* Aktualisierung Fehlerterm */
	 err -= es; 
	 if(err<0) {
		/* Fehlerterm wieder positiv (>=0) machen */
		err += el;
		/* Schritt in langsame Richtung, Diagonalschritt */
		x += ddx;
		y += ddy;
	 } else  {
		/* Schritt in schnelle Richtung, Parallelschritt */
		x += pdx;
		y += pdy;
	 }
	 array_INDEX2(out,int,0,t)=x;
	 array_INDEX2(out,int,1,t)=y;
  }

  /*------------------ /computation -------------------------*/

  return out;
}
