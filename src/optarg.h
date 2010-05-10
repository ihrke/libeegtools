/***************************************************************************
 *   Copyright (C) 2008/2009 by Matthias Ihrke                                  *
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

/**\file optarg.h
 \brief \ref status_stable Implementing optional argument lists.

 This feature can be used to conveniently pass arguments to 
 functions. This makes sense for generic functions that take
 one or more function pointers requiring different arguements each.

 The syntax is as follows:
 \code
 double *mydata = read_data();
 OptArgList *args = optarglist( "bottom_frequency=double,num_trials=int,data=double*",
                                0.5, 10, mydata );
 function_taking_optargs( args );
 free_optarglist( args ) // does not free any mydata 
 \endcode

 */
/*

 Every function that receives an optarglist is required to check for the 
 existence of "optarglist_free" within the parameterlist. If this parameter
 exists, the function is required to free the memory pointed to the
 argument. This enables the use of optarglists like this:
 \code
 function_taking_optargs( optarglist( "bottom_frequency=double,
                                       num_trials=int,
													data=double*,
													optarglist_free=int",
													0.5, 10, mydata, 1) );
 \endcode
 where the memory deallocation is handled by the called function.
 */
#ifndef OPTARG_H
# define OPTARG_H

#include "definitions.h"
#include <stdarg.h>

#define NO_OPTARGS NULL

#ifdef __cplusplus
extern "C" {
#endif
  /**\ingroup list
	*\{
	*/ 


  /** \brief Check and assign a scalar from OptArgList.
		Suppose you got optargs in a function and want to 
		check for variable "seed" of type "unsigned long",
		\code
		double x; // tmp 
		unsigned long seed;
		optarg_PARSE_SCALAR( optargs, "seed", seed, unsigned long, x );
		\endcode
	*/
#define optarg_PARSE_SCALAR( opts, label, var, type, tmp )	\
  if( optarglist_has_key( (opts), (label) ) ){					\
	 (tmp) = optarglist_scalar_by_key( (opts), (label) );		\
	 if( !isnan( (tmp) ) ) var=(type)(tmp);						\
  }


  /* ------------ List -------------- */
  struct list {
    void *content;
    struct list *next;
  };
  typedef struct list List;

  List* list_append( List *l );
  void  list_print( List *l, void(*print_content)(void*) );
  int   list_length( List *l );

  /**\ingroup optarg
	*\{
	*/ 
  OptArgList* optarglist( char *format, ... );

  bool        optarglist_has_key( OptArgList *list, const char *key );

  double      optarglist_scalar_by_key( OptArgList *list, const char *key );
  void*       optarglist_ptr_by_key   ( OptArgList *list, const char *key );
  OptArg*     optarglist_optarg_by_key( OptArgList *list, const char *key );

  OptArgList* optarglist_append_arg   ( OptArgList *list, OptArg *arg );

  void        optarglist_print( OptArgList *list, FILE *out );
  void        optarglist_free( OptArgList *list );

  OptArg*     optarg_scalar( const char *key, double scalar );
  OptArg*     optarg_ptr   ( const char *key, void *ptr );


  /*\}*/

#ifdef __cplusplus
}
#endif


#endif
