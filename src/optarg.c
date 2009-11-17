#include "optarg.h"
#include "helper.h"
#include <string.h>
#include <math.h>

List* list_append( List *l ){
  if( l==NULL ){
	 l = (List*)malloc( sizeof(List) );
	 return l;
  } else {
	 while( l->next!=NULL ) l=l->next;
	 l->next = (List*)malloc(sizeof(List));
	 return l->next;
  }
}

void  list_print( List *l, void(*print_content)(void*) ){
  if( l==NULL ){
	 printf("*\n");
	 return;
  }
  (*print_content)(l->content);
  printf( "->" );
  list_print( l->next, print_content );
}

int   list_length( List *l ){
  if( l==NULL )
	 return 0;
  else
	 return 1+list_length(l->next);
}

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

/** Creates an OptArgList out of one or many arguments

	 The Format string has the following convention:
	 <pre>
	 key1=double*,key2=void*,key3=int,key=double,...
	 </pre>

	 If you have at least one asterisk (*) in the specification, then
	 the corresponding argument needs to be a pointer. 
	 Internally, all pointer types are \c void*  and stored along
	 with the given specification.

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
OptArgList* optarglist( const char *format, ... ){
  va_list ap;
  char *tmp, *key, *type;
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

  L = (OptArgList*)malloc( sizeof( OptArgList ) );
  L->nargs=nargs;
  L->args=(OptArg*)malloc( nargs*sizeof( OptArg ) );
  dprintf("Found %i args\n", nargs );

  /* parse key=type information from format */
  fmt = (char*)malloc( strlen( format )*sizeof(char) );
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



  va_start (ap, format );         /* Initialize the argument list. */

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
		L->args[i].scalar=FALSE;
		L->args[i].data_ptr = va_arg (ap, void * );  
	 }
  }

  va_end (ap);                  /* Clean up. */

  free( fmt );
  return L;
}

/** does not free any data_ptr.
 */
void        optarglist_free( OptArgList *list ){
  free( list );
}


/** Return the value of the scalar argument or NaN if not available.
	 Check returned value with isnan()!
 */
double      optarglist_scalar_by_key( OptArgList *list, const char *key ){
  OptArg *arg;
  arg = optarglist_optarg_by_key( list, key );
  if( !arg ){
	 return NAN;
  }
  if( !arg->scalar ){
	 errprintf("You requested a scalar from a pointer argument\nYou get NaN instead :-)\n");
	 return NAN; 
  }
  return arg->data_scalar;
}

/** Return a pointer to the data_ptr field of the OptArg struct with the key from the
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

/** Return a pointer to an OptArg struct with the key from the
	 list. NULL if no such key has been found.
 */
OptArg*     optarglist_optarg_by_key( OptArgList *list, const char *key ){
  int i;
  for( i=0; i<list->nargs; i++ ){
	 if( !strcmp( list->args[i].key, key ) )
		return &(list->args[i]);
  }
  return NULL;
}
