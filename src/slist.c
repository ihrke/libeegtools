/***************************************************************************
 *   Copyright (C) 2008-2010 by Matthias Ihrke   *
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

#include "helper.h"
#include <string.h>
#include <stdlib.h>
#include "slist.h"

/** \brief initialize SingleList element.
 */
SingleList* slist_init( void *content ){
  SingleList *el;
  MALLOC( el, 1, SingleList );
  el->content=content;
  el->next=NULL;
  return el;
}

/** \brief append an element to the list.
 */
SingleList* slist_append( SingleList *l ){
  if( l==NULL ){
	 l = (SingleList*)malloc( sizeof(SingleList) );
	 return l;
  } else {
	 while( l->next!=NULL ) l=l->next;
	 l->next = slist_init(NULL);
	 return l->next;
  }
}
/** \brief return the idx'th element in the list.
	 
	 If it does not exist, NULL is returned.
	 \param l the head of the list
	 \param idx the index
	 \return the idx'th list node or NULL
*/
SingleList* slist_index( const SingleList *l, int idx ){
  SingleList *idxn=l;
  int i;
  for( i=0; i<idx; i++ ){
	 if( idxn==NULL ){
		return idxn;
	 }
	 idxn=idxn->next;
  }
  return idxn;
}

/** \brief print a list.
 */
void  slist_print( FILE *out, const SingleList *l, void(*print_content)(FILE*,void*) ){
  if( l==NULL ){
	 fprintf( out, "*\n" );
	 return;
  }
  (*print_content)(out,l->content);
  fprintf( out, "->" );
  slist_print( out, l->next, print_content );
}

/** \brief return the number of elements in a list.
 */
int   slist_length( const SingleList *l ){
  if( l==NULL )
	 return 0;
  else
	 return 1+slist_length(l->next);
}

/** \brief recursively free list.
	 The list->content pointer is NOT freed!
 */
void slist_free( SingleList *l ){
  if( l==NULL ) return;
  slist_free( l->next );
  free( l );
}
