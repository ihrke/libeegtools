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

 \section tempoptargs Temporary Optional Argument lists

 Every function that receives an optarglist is required to check for the 
 existence of "optarglist_free" within the parameterlist. If this parameter
 exists, the function is required to free the memory pointed to the
 argument. Also, the parameter is removed from the list and possibly 
 passed to other functions.

 This enables the use of optarglists like this:
 \code
 function_taking_optargs( optarglist( "bottom_frequency=double,
                                       num_trials=int,
													data=double*,
													optarglist_free=int",
													0.5, 10, mydata, 1) );
 \endcode
 or even shorter like
 \code
 function_taking_optargs( optarglisttmp( "bottom_frequency=double,
                                        num_trials=int,
													 data=double*",												
												 	0.5, 10, mydata) );
 \endcode
 where the memory deallocation is handled by the called function.

 \par For writing a function that takes an OptArgList* as an argument.

	 The function optarlist_remove_freeflag() 
	 removes the optarglist_free-flag from the list
	 and returns a new list without the flag. The old list is free'd 
	 in case that the optarglist_free argument is found. 

	 The removedflag-flag is set to FALSE if it is the same list,
	 to TRUE if it is a new (truncated) list.

	 Usage is as follows:
	 \code
	 void testfct( OptArgList *opts ){	 
	    // overwrite the pointer
		 bool freeflag_removed=FALSE;
	    opts=optarglist_remove_freeflag( opts, &freeflag_removed );

		 // do some stuff including passing the argument to some other function
		 call_other_function( opts );
		 
		 // clean up in case the flag was indeed removed
		 if( freeflag_removed )
		    optarglist_free( opts );
	 }
	 \endcode

*/
#ifndef OPTARG_H
# define OPTARG_H

#include "definitions.h"
#include "slist.h"
#include <stdarg.h>

#define NO_OPTARGS NULL

#ifdef __cplusplus
extern "C" {
#endif
  /*-----------------------------------------------------------
	 - Optional Arguments
	 ---------------------------------------------------------*/
  /**\brief An optional Argument. */
  typedef struct {
	 char    key[MAX_LABEL_LENGTH]; /**<\brief The key to refer to the argument */
	 bool    scalar;                /**<\brief flag whether arg is scalar or pointer*/
	 char    type[MAX_LABEL_LENGTH];/**<\brief e.g. double* ... */
	 void    *data_ptr;             /**<\brief data if scalar=FALSE */
	 double  data_scalar;           /**<\brief data if scalar=TRUE */
  } OptArg;

  /**\brief A list of optional Arguments. */
  typedef struct {
	 int nargs; /**<\brief number of arguments in the list */
	 OptArg *args; /**<\brief memory for the arguments */
  } OptArgList;

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
	 if( !isnan( (tmp) ) ) (var)=(type)(tmp);						\
  }

 /** \brief Check and assign a void-ptr from OptArgList.
		Suppose you got optargs in a function and want to 
		check for variable "distfct" of type "DistFct",
		\code
		void *ptr; // tmp 
		DistFct distf;
		optarg_PARSE_PTR( optargs, "distfct", distf, DistFct, ptr );
		\endcode
	*/
#define optarg_PARSE_PTR( opts, label, var, type, tmp )		\
  if( optarglist_has_key( (opts), (label) ) ){					\
	 (tmp) = optarglist_ptr_by_key( (opts), (label) );			\
	 if( (tmp) ) (var)=(type)(tmp);									\
  }

  OptArgList* optarglist( char *format,  ... );
  OptArgList* optarglisttmp( char *format, ... );

  bool        optarglist_has_key( OptArgList *list, const char *key );

  double      optarglist_scalar_by_key( OptArgList *list, const char *key );
  void*       optarglist_ptr_by_key   ( OptArgList *list, const char *key );
  OptArg*     optarglist_optarg_by_key( OptArgList *list, const char *key );

  OptArgList* optarglist_append_arg   ( OptArgList *list, OptArg *arg );
  OptArgList* optarglist_delete_arg   ( OptArgList *list, OptArg *arg );

  void        optarglist_print( OptArgList *list, FILE *out );
  void        optarglist_free( OptArgList *list );

  OptArgList* optarglist_remove_freeflag( OptArgList *list, bool *removedflag );

  OptArg*     optarg_scalar( const char *key, double scalar );
  OptArg*     optarg_ptr   ( const char *key, void *ptr );

#ifdef __cplusplus
}
#endif


#endif
