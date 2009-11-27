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

#include<stdlib.h>
#include<check.h>
#include<math.h>
#include "check_all.h"
#include "helper.h"

     
START_TEST (test_strip_blank)
{
  char s1[40]=" thi s is a blank st ri n g ";
  char s2[40]="  more blanks     at the end     ";

  string_strip_blanks( (char*)s1 );
  string_strip_blanks( (char*)s2 );

  fail_if( strcmp( s1, "thisisablankstring" ), "FAILURE: got '%s'\n", s1 );
  fail_if( strcmp( s2, "moreblanksattheend" ), "FAILURE: got '%s'\n", s2 );

}
END_TEST



Suite * init_other_suite (void){
  Suite *s = suite_create ("Other/Misc-Functions");

  TCase *tc_core = tcase_create ("OtherCore");
  tcase_add_test (tc_core, test_strip_blank );

  tcase_set_timeout(tc_core, 20);
  suite_add_tcase (s, tc_core);
  return s;
}

