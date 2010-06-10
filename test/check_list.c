/*
** check_list.c
** 
** Unit-Test suite using libcheck
**
**
** Made by (Matthias Ihrke)
** Login   <mihrke@localhost>
** 
*/

#include <stdlib.h>
#include <check.h>
#include "check_all.h"
#include <math.h>
#include "definitions.h"
#include "helper.h"
#include "slist.h"

void print_int( FILE *out,void *c){
  int co=*(int*)c;
  fprintf(out,"%i",co);
}

START_TEST (test_append)
{
  int c[10]={1,2,3,4,5,6,7,8,9,10};
  SingleList *l=slist_init( (void*)c );
  SingleList *tmp;
  int i;
  for( i=1; i<10; i++ ){
	 tmp=slist_append( l );
	 tmp->content=(void*)(c+i);
  }

  tmp = l;
  i=1;
  while( tmp->next!=NULL ){
	 fail_if( *(int*)tmp->content != i++ );
	 tmp=tmp->next;
  }

  /* slist_print( stderr, l, print_int ); */
  slist_free( l );
}
END_TEST

START_TEST (test_index)
{
  int c[10]={1,2,3,4,5,6,7,8,9,10};
  SingleList *l=slist_init( (void*)c );
  SingleList *tmp;
  int i;
  for( i=1; i<10; i++ ){
	 tmp=slist_append( l );
	 tmp->content=(void*)(c+i);
  }

  for( i=0; i<10; i++ ){
	 tmp=slist_index( l, i );
	 fail_if( *(int*)tmp->content != (i+1) );
  }

  /* slist_print( stderr, l, print_int ); */
  slist_free( l );
}
END_TEST

START_TEST (test_length)
{
  int c[10]={1,2,3,4,5,6,7,8,9,10};
  SingleList *l=slist_init( (void*)c );
  SingleList *tmp;
  int i;
  for( i=1; i<10; i++ ){
	 fail_unless( slist_length(l)==i );
	 tmp=slist_append( l );
	 tmp->content=(void*)(c+i);
	 
  }
  /* slist_print( stderr, l, print_int ); */
  slist_free( l );
}
END_TEST

/* template
START_TEST (test_)
{
}
END_TEST
*/

Suite * init_list_suite (void){
  Suite *s = suite_create ("List-Functions");

  TCase *tc_core = tcase_create ("ListCore");
  tcase_add_test (tc_core, test_append );
  tcase_add_test (tc_core, test_index );
  tcase_add_test (tc_core, test_length );

  tcase_set_timeout(tc_core, 20);
  suite_add_tcase (s, tc_core);
  return s;
}

