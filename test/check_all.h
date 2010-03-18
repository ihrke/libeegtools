
/*
** check_all.h
** 
** Made by Matthias Ihrke
** Login   <mihrke@localhost>
** 
** Started on  Thu Oct 11 13:26:09 2007 Matthias Ihrke
** Last update Thu Oct 11 13:40:35 2007 Matthias Ihrke
*/

#ifndef   	CHECK_ALL_H_
# define   	CHECK_ALL_H_
#include<check.h>
#include<stdlib.h>
#include<math.h>
#include<stdio.h>

Suite * init_denoising_suite (void);
Suite * init_checkfct_suite  (void);
Suite * init_distance_suite  (void);
Suite * init_hmm_suite       (void);
Suite * init_list_suite      (void);
Suite * init_other_suite     (void);
Suite * init_nonlinear_suite (void);
Suite * init_array_suite     (void);

#endif 	    /* !CHECK_ALL_H_ */
