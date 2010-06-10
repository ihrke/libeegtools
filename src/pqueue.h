/***************************************************************************
 *   Copyright (C) 2010 by Matthias Ihrke                                  *
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
 *   aint with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/**\file pqueue.h
 \brief \ref status_stable Priority queue.

Smaller priority values indicate higher priority.

The implementation is based on a single-linked list. 
Popping elements is therefore in O(1), inserting in O(N).

Usage is as follows:

\code
  PriorityQueue *pq=pq_init();
  int stuff[] = {1, 2, 3, 4};
  
  pq_insert( pq, (void*)&(stuff[0]), 10 );
  ...
  pq_insert( pq, (void*)&(stuff[3]), 1 );

  int priorel=*((int*)pq_pop( pq ));
  int nextprior=*((int*)pq_pop( pq ));
  ...

  pq_free( pq );
\endcode
	
 */
#ifndef PQUEUE_H
# define PQUEUE_H
#include "definitions.h"
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

  /*-----------------------------------------------------------
	 - Priority Queue
	 ---------------------------------------------------------*/

  /** \brief a node in the PQueue.
	*/
  struct pqueue {
	 double priority;
	 void *content;
	 struct pqueue *next;
  };
  typedef struct pqueue PQnode;

  /** \brief the head of the PQueue.
	*/
  typedef struct {
	 PQnode *root;
  } PriorityQueue;


  /* -------------- FUNCTIONS ---------------- */
  PriorityQueue* pq_init();
  void    pq_insert( PriorityQueue *pq, void *c, double prior );
  void*   pq_pop( PriorityQueue *pq );
  PQnode* pq_pop_node( PriorityQueue *pq );
  void    pq_print( FILE *out, PriorityQueue *pq );

  void    pq_free( PriorityQueue *pq );

  /**\cond PRIVATE */
  PQnode* pqnode_init( void *c, double priority );
  void    pqnode_print( FILE *out, PQnode *N );
  void    pqnode_free( PQnode *n );
  /**\endcond */


#ifdef __cplusplus
}
#endif

#endif /* PQUEUE_H */
