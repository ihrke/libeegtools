/***************************************************************************
 *   Copyright (C) 2008-2010 by Matthias Ihrke                                  *
 *   mihrke@uni-goettingen.de                                              *
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

/**\file slist.h
 \brief \ref status_stable Implementing a single-linked list.

 */
#ifndef SLIST_H
#define SLIST_H

#include <stdio.h>
#include "definitions.h"

#ifdef __cplusplus
extern "C" {
#endif

  struct slist {
    void *content;    /**< the content of this list's entry */
    struct list *next; /**< pointer to next list entry */
  }; 

  /** \brief a single-linked list.
	*/
  typedef struct slist SingleList;

  SingleList* slist_init  ( void *content );
  SingleList* slist_append( SingleList *l );
  void        slist_print ( FILE *out, const SingleList *l, void(*print_content)(FILE*,void*) );
  int         slist_length( const SingleList *l );
  SingleList* slist_index ( const SingleList *l, int idx );
  void        slist_free  ( SingleList *l );

#ifdef __cplusplus
}
#endif


#endif /* SLIST_H */
