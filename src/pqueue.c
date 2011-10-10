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

#include <stdlib.h>

#include "pqueue.h"
#include "helper.h"

/** TEMPLATE
	 \param
	 \param
	 \return
*/
PriorityQueue* pq_init( ) {
  PriorityQueue *pq=(PriorityQueue*)malloc(sizeof(PriorityQueue));
  pq->root=NULL;
  return pq;
}

PQnode* pqnode_init( void *c, double priority ){
  PQnode *pq;
  pq = (PQnode*)malloc(sizeof(PQnode));
  pq->priority=priority;
  pq->content = c;
  pq->next = NULL;
  return pq;
}

void  pq_insert( PriorityQueue *pq, void *c, double prior ){
  PQnode *node=pqnode_init( c, prior );
  PQnode *it=pq->root;
  
  dprintf("Inserting node %p with priority %f\n", c, prior );
  if( it==NULL || prior<it->priority ){
	 node->next=pq->root;
	 pq->root=node;
	 return;
  }
  
  while( 1 ){
	 if( it->next==NULL ){
		it->next = node;
		return;
	 } else if( prior<it->next->priority ){
		node->next = it->next;
		it->next = node;
		return;
	 } else {
		it=it->next;
	 }
  }
}

/** pop a node in the PQ. use it if you need access to
	 the priority value as well. You will need to free
	 the node after usage yourself;
	 \param pq the PQ
*/
PQnode* pq_pop_node( PriorityQueue *pq ){
  PQnode *tmp;
  if( pq->root==NULL ){
	 warnprintf("Popping from empty PQ\n");
	 return NULL;
  }
  tmp = pq->root;
  pq->root = tmp->next;
  tmp->next=NULL;
  return tmp;
}

/** pop the data with highest priority in the PQ.
	 if you need access to the priority as well, use pq_pop_node().
	 The corresponding node is removed from the PQ;
	 \param pq the PQ
*/
void* pq_pop( PriorityQueue *pq ){
  PQnode *tmp;
  void *content;

  tmp = pq_pop_node( pq );
  content = tmp->content;
  pqnode_free( tmp );

  return content;
}


void pq_print( FILE *out, PriorityQueue *pq ){
  fprintf( out, "PQ(%p)=[ ", pq );
  pqnode_print( out, pq->root );
  fprintf( out, " ]\n");
}
/**\cond PRIVATE */
void pqnode_print( FILE *out, PQnode *N ){
  if( N==NULL ){
	 return;
  }
  fprintf( out, "(%p,%.2f) ", N->content, N->priority );
  pqnode_print( out, N->next );
}

void pqnode_free( PQnode *n ){
  if( !n ) return;
  if( n->next ){
	 pqnode_free( n->next );
  }
  free( n );
}
/**\endcond */

void  pq_free( PriorityQueue *pq ){
  pqnode_free( pq->root );
  free( pq );
  return;
}
