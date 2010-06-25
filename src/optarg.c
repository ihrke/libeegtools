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
#include "optarg.h"
#include "helper.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>

/** \brief print an optarglist.
 */
void optarglist_print( OptArgList *list, FILE *out ){
  int i;
  fprintf( out, "OptArgList '%p' with %i items:\n", list, list->nargs );
  for( i=0; i<list->nargs; i++ ){
	 fprintf( out, " %i: %s(%s) is %s: ", i+1, list->args[i].key, list->args[i].type,
				 list->args[i].scalar?"scalar":"pointer");
	 if( list->args[i].scalar ){
		fprintf( out, "%f\n", list->args[i].data_scalar );
	 } else {
		fprintf( out, "%p\n", list->args[i].data_ptr );
	 }
  }
}

/**\cond PRIVATE
 */

OptArgList* optarglist_hidden( char *format, va_list ap ){
  char *key, *type;
  char *tmp;
  char *fmt;
  int nargs, i;
  OptArgList *L;

  nargs=0;
  tmp=format;

  /* count '=' */
  while( (tmp=strchr( tmp, '=' ))!=NULL ){
	 tmp++;
	 nargs++;
  }

  if( sizeof(char*)!=sizeof(void*) ||
		sizeof(double*)!=sizeof(void*) ){
	 errprintf("System Error! For this function to work, all pointers must be of the same size!\n");
	 return NULL;
  }

  L = (OptArgList*)malloc( sizeof( OptArgList ) );
  L->nargs=nargs;
  L->args=(OptArg*)malloc( nargs*sizeof( OptArg ) );
  dprintf("Found %i args\n", nargs );

  /* parse key=type information from format */
  fmt = (char*)malloc( (strlen( format )+1)*sizeof(char) );
  strcpy( fmt, format );
  key=fmt;
  type=fmt;
  for( i=0; i<nargs; i++ ){
	 if( (tmp=strchr( key, '=' )) ){
		*tmp='\0';
	 }
	 dprintf("key=%s\n", key);
	 strcpy( L->args[i].key, key );
	 string_strip_blanks( L->args[i].key );

	 type=tmp+1;
	 if( (tmp=strchr( type, ',' ) ) ){
		*tmp='\0';
		key = tmp+1;
	 }

	 dprintf("type=%s\n", type );
	 strcpy( L->args[i].type, type );
	 string_strip_blanks( L->args[i].type );
  }
  dprintf("Parsing done\n");

  for (i = 0; i < nargs; i++){
	 if( !strchr( L->args[i].type, '*' ) ){ /* it's scalar */
		L->args[i].scalar=TRUE;
		if( !strcmp( L->args[i].type, "char" ) ){
		  L->args[i].data_scalar = (double) va_arg (ap, int );  
		} else if( !strcmp( L->args[i].type, "short" ) ){
		  L->args[i].data_scalar = (double) va_arg (ap, int );  
		} else if( !strcmp( L->args[i].type, "int" ) ){
		  L->args[i].data_scalar = (double)va_arg (ap, int );  
		} else if( !strcmp( L->args[i].type, "long" ) ){
		  L->args[i].data_scalar = (double) va_arg (ap, long );  
		} else if( !strcmp( L->args[i].type, "float" ) ){
		  L->args[i].data_scalar = (double) va_arg (ap, double );  
		} else if( !strcmp( L->args[i].type, "double" ) ){
		  L->args[i].data_scalar = (double) va_arg (ap, double );  
		} else {
		  errprintf("Sorry, do not know scalar data type '%s'\n"
						"This function is probably going to break because of this.\n"
						"Please correct the error and restart.\n",
						L->args[i].type );
		  break;
		}
	 } else {
		dprintf("There is a pointer in the %i'th argument\n", i+1 );
		L->args[i].scalar=FALSE;
		L->args[i].data_ptr = (void*)va_arg(ap, void *);  
		dprintf("got a pointer: %p\n", L->args[i].data_ptr );
	 }
  }


  free( fmt );
  return L;
}

/**\endcond */

/**\brief Creates an OptArgList including a optarglist_free arg.

	 The Format string has the following convention:
	 <pre>
	 key1=double*,key2=void*,key3=int,key=double,...
	 </pre>

	 If you have at least one asterisk (*) in the specification, then
	 the corresponding argument needs to be a pointer. 
	 Internally, all pointer types are \c void*  and stored along
	 with the given specification.

	 Be VERY careful to typecast all variables to match the type you indicated.
	 For example, if you have a float-variable and pass it as double, you need to 
	 typecast to (double).

	 Only "basic" scalar values (without asterisk) are
	 supported:
	 <ul>
	 <li> char
	 <li> short
	 <li> int
	 <li> long
	 <li> float
	 <li> double
	 </ul>

	 Blanks/Newlines etc are stripped from key and type, so
	 <pre>
	 key = double , key2=void*
	 </pre>
	 are ok.
*/
OptArgList* optarglisttmp( char *format, ... ){
  va_list ap;
  va_start (ap, format );         /* Initialize the argument list. */
  OptArgList *o=optarglist_hidden( format, ap );
  va_end(ap);
  
  OptArg f;
  sprintf(f.key, "optarglist_free");
  f.scalar=TRUE;
  sprintf(f.type, "int");
  f.data_ptr=NULL;
  f.data_scalar=1.0;
  OptArgList *no=optarglist_append_arg( o, &f );
  optarglist_free( o );
  
  return no;
}

/** \brief Creates an OptArgList out of one or many arguments.

	 The Format string has the following convention:
	 <pre>
	 key1=double*,key2=void*,key3=int,key=double,...
	 </pre>

	 If you have at least one asterisk (*) in the specification, then
	 the corresponding argument needs to be a pointer. 
	 Internally, all pointer types are \c void*  and stored along
	 with the given specification.

	 Be VERY careful to typecast all variables to match the type you indicated.
	 For example, if you have a float-variable and pass it as double, you need to 
	 typecast to (double).

	 Only "basic" scalar values (without asterisk) are
	 supported:
	 <ul>
	 <li> char
	 <li> short
	 <li> int
	 <li> long
	 <li> float
	 <li> double
	 </ul>

	 Blanks/Newlines etc are stripped from key and type, so
	 <pre>
	 key = double , key2=void*
	 </pre>
	 are ok.
 */
OptArgList* optarglist( char *format, ... ){
  va_list ap;
  va_start (ap, format );         /* Initialize the argument list. */
  OptArgList *o=optarglist_hidden( format, ap );
  va_end(ap);

  return o;
}

/** \brief free an optarglist.

	 does not free any data_ptr.
 */
void        optarglist_free( OptArgList *list ){
  free( list->args );
  free( list );
}


/** \brief Return the value of the scalar argument.

	 or NaN if not available.
	 Check returned value with isnan()!
 */
double      optarglist_scalar_by_key( OptArgList *list, const char *key ){
  OptArg *arg;
  arg = optarglist_optarg_by_key( list, key );
  if( !arg ){
	 errprintf("There is no key: '%s'\n", key );
	 return NAN;
  }
  if( !arg->scalar ){
	 errprintf("You requested a scalar from a pointer argument\nYou get NaN instead :-)\n");
	 return NAN; 
  }
  dprintf("arg->data_scalar=%f\n", arg->data_scalar );
  return arg->data_scalar;
}

/** \brief Return a pointer to the data_ptr field of the OptArg struct.

	 with the key from the
	 list. NULL if no such key has been found.
 */
void*       optarglist_ptr_by_key   ( OptArgList *list, const char *key ){
  OptArg *arg;
  arg = optarglist_optarg_by_key( list, key );
  if( !arg ) return NULL;
  if( arg->scalar ){
	 errprintf("You requested a data-pointer from a scalar argument\nYou get NULL instead :-)\n");
	 return NULL; 
  }
  return arg->data_ptr;
}

/** \brief Return a pointer to an OptArg struct with the key from the
	 list. 

	 NULL if no such key has been found.
 */
OptArg*     optarglist_optarg_by_key( OptArgList *list, const char *key ){
  int i;
  if( !list ) return NULL;
  for( i=0; i<list->nargs; i++ ){
	 if( !strcmp( list->args[i].key, key ) )
		return &(list->args[i]);
  }
  return NULL;
}

/** \brief return TRUE if key is found in list, else FALSE.
 */
bool     optarglist_has_key( OptArgList *list, const char *key ){
  int i;
  if( !list ) return FALSE;
  for( i=0; i<list->nargs; i++ ){
	 if( !strcmp( list->args[i].key, key ) )
		return TRUE;
  }
  return FALSE;
}

/** \brief create a new optarglist that is the former optarglist with the new arg appended.

	 The caller is responsible for freeing the old list.

	 \param list the pointer to the Optarglist
	 \return the new list (the caller is responsible for freeing the old one)
 */
OptArgList*     optarglist_append_arg   ( OptArgList *list, OptArg *arg ){
  OptArgList *new;
  
  new = (OptArgList*)malloc( sizeof( OptArgList ) );
  new->nargs=list->nargs+1;
  new->args=(OptArg*)malloc( new->nargs*sizeof( OptArg ) );
  memcpy( new->args, list->args, list->nargs * sizeof( OptArg ) );
  memcpy( &(new->args[new->nargs-1]), arg, sizeof( OptArg ) );
  return new;
}

/**\brief delete an argument from the argument list.

	The function returns a new OptArgList containing 
	which is the old list minus the argument. The arguments 
	copied!

	The caller is responsible for freeing the old list. 

	\param list the list
	\param arg the argument to delete
	\return the new list (the caller is responsible for freeing the old one)
 */
OptArgList* optarglist_delete_arg   ( OptArgList *list, OptArg *arg ){
  int idx=-1;
  int i,j;
  for( i=0; i<list->nargs; i++ ){
	 if( &(list->args[i])==arg ){
	 	idx=i;
	 }
  }
  if( idx<0 ){
	 warnprintf("list did not contain the optional argument\n");
	 return;
  }
  
  OptArgList *L = (OptArgList*)malloc( sizeof( OptArgList ) );
  L->nargs=list->nargs-1;
  L->args=(OptArg*)malloc( L->nargs*sizeof( OptArg ) );
  for( i=0,j=0; i<list->nargs; i++ ){
	 if( i==idx ) continue;
	 memcpy( &(L->args[i]), &(list->args[i]), sizeof(OptArg) );
	 j++;
  }

  return L;
}
/** \brief remove the optarglist_free-flag from the list.
	 
	 This function is only interesting for you, if you want to write an
	 own function that takes an OptArgList* as an argument.

	 The function removes the optarglist_free-flag from the list
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
	 \endcode

	 \param list the list
	 \param removedflag is set to FALSE if it is the same list,
	         to TRUE if it is a new (truncated) list.
	 \return a pointer to the new (or the old) list
 */
OptArgList* optarglist_remove_freeflag( OptArgList *list, bool *removedflag ){
  OptArg *freearg;
  OptArgList *nlist=list;
  *removedflag=FALSE;
  if( optarglist_has_key( list, "optarglist_free" ) ){
	 freearg=optarglist_optarg_by_key( list, "optarglist_free" );
	 nlist=optarglist_delete_arg( list, freearg );
	 optarglist_free( list );
	 *removedflag=TRUE;
  } 
  
  return nlist;
}

/** \brief create a new OptArg* containing a scalar value.
 */
OptArg*     optarg_scalar( const char *key, double scalar ){
  OptArg *p = (OptArg*)malloc(sizeof(OptArg) );
  strncpy( p->key, key, MAX_LABEL_LENGTH );
  p->scalar = TRUE;
  p->data_scalar = scalar;
  sprintf( p->type, "double" );
  p->data_ptr = NULL;
  return p;
}

/** \brief create a new OptArg* containing a pointer.
 */
OptArg*     optarg_ptr   ( const char *key, void *ptr ){
  OptArg *p = (OptArg*)malloc(sizeof(OptArg) );
  strncpy( p->key, key, MAX_LABEL_LENGTH );
  p->scalar = FALSE;
  p->data_ptr = ptr;
  sprintf( p->type, "void*" );
  p->data_scalar = 0.0/0.0;
  return p;
}
