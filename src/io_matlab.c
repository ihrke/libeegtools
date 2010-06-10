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

#include "io_matlab.h"

#ifdef MATIO
#include <matio.h>

/* ------------------------------- READER ----------------------------- */

double get_double_from_struct_field( matvar_t *eeg, const char *name, int struct_array_index );

/** \brief read EEG struct from EEGlab (MATLAB) .set-file.

	 This reader uses MatIO to parse EEGlab files and construct
	 a LibEEGTools EEG struct. It was developed with EEGlab version 6.01b.
	 
	 Currently, the EEGlab set must be "epoched", that means that the data must be
	 three-dimensional with channels x n x trials. Continuous data is not yet supported
	 and the reader will fail.
	 
	 The reader currently works only if the storage is in single .set file (as opposed to a 
	 pair of .set and .dat file).
	 
	 \param file eeglab .set file
	 \return the EEG struct
*/
EEG* read_eeglab_file( const char *file ){
  mat_t *mat;
  matvar_t *meeg; /* contains the struct 'EEG' from EEGlab */
  matvar_t *tmp, *tmp2, *event, *epoch;  
  int nfields; 
  int c,i,j;
  EEG *eeg;
  int nbchan, ntrials, n;

  dprintf("Reading file: '%s'\n", file);
  mat = Mat_Open( file, MAT_ACC_RDONLY);
  if( !mat ){
	 errprintf("Error opening MATLAB file '%s'\n", file );
	 return NULL;
  }
  meeg = Mat_VarRead( mat, "EEG" );
  if( meeg->class_type!=MAT_C_STRUCT ){
	 errprintf("EEG does not appear to be a struct\n" );
	 return NULL;
  }
  nfields = Mat_VarGetNumberOfFields( meeg );
#ifdef DEBUG
  dprintf( "There are %i fields in the EEG struct\n", nfields );
  for( i=1; i<=nfields; i++ ){ /* MATLAB 1-relative indices */
	 tmp = Mat_VarGetStructField( meeg, &i, BY_INDEX, 0 );
	 dprintf("Found field: '%s'\n", tmp->name);
  }
#endif
  
  /* dimensions */
  nbchan = (int)get_double_from_struct_field( meeg, "nbchan",0 );
  ntrials= (int)get_double_from_struct_field( meeg, "trials",0 );
  n      = (int)get_double_from_struct_field( meeg, "pnts",0 );
  dprintf("dim=(%i,%i,%i)\n", nbchan,ntrials,n);
  eeg = eeg_init( nbchan, ntrials, n );

  /* filename  */
  eeg->filename=(char*)malloc( (strlen(file)+2)*sizeof(char) );
  strcpy( eeg->filename, file );

  /* comments */
  tmp = Mat_VarGetStructField( meeg, "comments", BY_NAME, 0 );
  eeg->comment=(char*)malloc( tmp->dims[0]*sizeof(char) );
  for( i=0; i<tmp->dims[0]+1; i++ ){
	 eeg->comment[i]='\0';
  }
  memcpy( eeg->comment, tmp->data, tmp->dims[0]*sizeof(char));

  /* sampling rate */
  eeg->sampling_rate = get_double_from_struct_field( meeg, "srate", 0);

  /* times */
  tmp = Mat_VarGetStructField( meeg, "times", BY_NAME, 0 );
  if( tmp->dims[1]==0 && ntrials == 1){ // continuous data
	 dprintf("Continuous data, skipping times-array\n");
  } else if( tmp->dims[1]!=n ){
	 errprintf("times-array should be of length n: %i != %i\n", tmp->dims[1], n );
	 eeg_free( eeg );
	 return NULL;
  } else {
	 if( tmp->data_size != sizeof(double) ){
		errprintf("times is not double format, %i!=%li\n", tmp->data_size, sizeof(double));
		eeg_free( eeg );
		return NULL;	 
	 }
	 eeg->times=array_new2( DOUBLE, 1, n );
	 memcpy(eeg->times->data, tmp->data, n*sizeof(double) );
  }

  /* channel info */
  eeg->chaninfo = (ChannelInfo*)malloc( nbchan*sizeof(ChannelInfo) );
  tmp = Mat_VarGetStructField( meeg, "chanlocs", BY_NAME, 0 );
  dprintf("chanlocs: %i,%i\n", tmp->dims[0], tmp->dims[1]);
  
  for( i=0; i<nbchan; i++ ){
	 tmp2 = Mat_VarGetStructField( tmp, "labels", BY_NAME, i );
	 eeg->chaninfo[i].num = i;
	 eeg->chaninfo[i].num_chans = nbchan;
	 strcpy(eeg->chaninfo[i].label, (char*)tmp2->data);
	 eeg->chaninfo[i].x = get_double_from_struct_field( tmp, "X", i);
	 eeg->chaninfo[i].y = get_double_from_struct_field( tmp, "Y", i);
	 eeg->chaninfo[i].z = get_double_from_struct_field( tmp, "Z", i);
  }
  
  /* data */
  tmp = Mat_VarGetStructField( meeg, "data", BY_NAME, 0 );
  if( ntrials==1 ) { // continuous data 
	 if( tmp->dims[0]!=nbchan || tmp->dims[1]!=n ){
		errprintf("(nbchan,n)=(%i,%i), should be (%i,%i)\n",
					 tmp->dims[0], tmp->dims[1], nbchan, n );
		eeg_free( eeg );
		return NULL;
	 }
  } else  if( tmp->dims[0]!=nbchan || tmp->dims[1]!=n || tmp->dims[2]!=ntrials ){
	 errprintf("(nbchan,ntrials,n)=(%i,%i,%i), should be (%i,%i,%i)\n",
				  tmp->dims[0], tmp->dims[2], tmp->dims[1], nbchan, ntrials, n );
	 eeg_free( eeg );
	 return NULL;
  }
  if( tmp->data_size != sizeof(float) ){
	 errprintf("data is not in float format, sizeof(data)=%i, sizeof(float)=%li\n", 
				  tmp->data_size, sizeof(float));
	 eeg_free( eeg );
	 return NULL;	 
  } 
  float x;
  for( c=0; c<nbchan; c++ ){
	 for( i=0; i<ntrials; i++ ){
		for( j=0; j<n; j++ ){
		  x=((float*)tmp->data)[ c + (j*nbchan) + (i*n*nbchan) ];
		  array_INDEX3(eeg->data,double,c,i,j) = (double)x;
		}
	 }
  }
  
  /* TODO CONTINUE HERE */

  /* /\* markers *\/ */
  /* epoch = Mat_VarGetStructField( meeg, "epoch", BY_NAME, 0 ); */
  /* if( epoch->dims[0] == 0 || epoch->dims[1] < ntrials ){ */
  /* 	 warnprintf("no epoch field, or wrong dimensions (%i,%i), skipping...\n",  */
  /* 				  epoch->dims[0],epoch->dims[1]); */
  /* } else { */
  /* 	 eeg->nmarkers = (unsigned int*) malloc( ntrials*sizeof(unsigned int) ); */
  /* 	 eeg->markers = (unsigned int**) malloc( ntrials*sizeof(unsigned int*) ); */
  /* 	 eeg->marker_labels = (char***) malloc( ntrials*sizeof(char**) ); */
	 
  /* 	 for( i=0; i<ntrials; i++ ){ */
  /* 		tmp  = Mat_VarGetStructField( epoch, "eventlatency", BY_NAME, i ); */
  /* 		tmp2 = Mat_VarGetStructField( epoch, "eventtype", BY_NAME, i ); */
  /* 		dprintf("%i, %i\n", i, tmp->dims[1]); */
  /* 		eeg->nmarkers[i] = tmp->dims[1]; */
  /* 		eeg->markers[i] = (unsigned int*) malloc( eeg->nmarkers[i]*sizeof(unsigned int) ); */
  /* 		eeg->marker_labels[i] = (char**) malloc( eeg->nmarkers[i]*sizeof(char*) ); */
  /* 		for( j=0; j<eeg->nmarkers[i]; j++ ){ */
  /* 		  /\* label *\/ */
  /* 		  event = Mat_VarGetCell( tmp2, j ); /\* MATLAB index *\/ */
  /* 		  eeg->marker_labels[i][j] = (char*)malloc( (strlen((char*)event->data)+1)*sizeof(char) ); */
  /* 		  strcpy( eeg->marker_labels[i][j], (char*)event->data ); */
  /* 		  /\* latency *\/ */
  /* 		  event = Mat_VarGetCell( tmp, j ); /\* MATLAB index *\/ */
  /* 		  eeg->markers[i][j] = closest_index( eeg->times, n, ((double*)event->data)[0] ); */
  /* 		} */
  /* 	 } */
  /* } */
  dprintf("Finished reading '%s'\n", file );

  return eeg;
}

/** \brief read a MATLAB-array from a .mat file.
	 
	 \todo Currently, the returned array is always DOUBLE.

	 \param file the .mat file
	 \param varname the name of the variable in the .mat file; can be NULL, in this case the
	       first variable is read.
	 \return the Array-struct or NULL in case an error occurred
 */
Array* read_array_matlab( const char *file, const char *varname ){ 
  mat_t *mfile;
  matvar_t *marr=NULL;

  dprintf("Reading variable '%s' from file: '%s'\n", varname, file);
  mfile = Mat_Open( file, MAT_ACC_RDONLY);
  if( !mfile ){
	 errprintf("Error opening MATLAB file '%s'\n", file );
	 return NULL;
  }
  if( varname ){
	 marr = Mat_VarRead( mfile, varname );
  } else {
	 marr = Mat_VarReadNext( mfile );
  }
  Mat_Close( mfile );
  if( !marr ){
	 errprintf("Something is wrong, could not read variable\n");
	 return NULL;
  }
  Array *out=array_new( DOUBLE, marr->rank, marr->dims );

  ulong i; 
  uint *index=(uint*)malloc( out->ndim*sizeof(uint));
  for( i=0; i<marr->nbytes/marr->data_size; i++ ){
	 array_calc_colindex( i, out->size, out->ndim, index );
	 
	 switch( marr->data_type ){
	 case MAT_T_INT8:
		*((double*)array_index(out,index))=(double)(*(int8_t*)(marr->data+(i*marr->data_size))); 
		break;
	 case MAT_T_UINT8:
		*((double*)array_index(out,index))=(double)(*(uint8_t*)(marr->data+(i*marr->data_size)));
		break;
	 case MAT_T_INT16:
		*((double*)array_index(out,index))=(double)(*(int16_t*)(marr->data+(i*marr->data_size)));
		break;
	 case MAT_T_UINT16:
		*((double*)array_index(out,index))=(double)(*(uint16_t*)(marr->data+(i*marr->data_size)));
		break;
	 case MAT_T_INT32:
		*((double*)array_index(out,index))=(double)(*(int32_t*)(marr->data+(i*marr->data_size)));
		break;
	 case MAT_T_UINT32:
		*((double*)array_index(out,index))=(double)(*(uint32_t*)(marr->data+(i*marr->data_size)));
		break;
	 case MAT_T_INT64:
		*((double*)array_index(out,index))=(double)(*(int64_t*)(marr->data+(i*marr->data_size)));
		break;
	 case MAT_T_UINT64:
		*((double*)array_index(out,index))=(double)(*(uint64_t*)(marr->data+(i*marr->data_size)));
		break;
	 case MAT_T_SINGLE:
		*((double*)array_index(out,index))=(double)(*(float*)(marr->data+(i*marr->data_size)));
		break;
	 case MAT_T_DOUBLE:
		*((double*)array_index(out,index))=(double)(*(double*)(marr->data+(i*marr->data_size)));
		break;
	 default:
		errprintf("Unknown Data-Type in MATLAB-file!\n");
		break;
	 }
  }
  free(index);
  
  return out;
}

double get_double_from_struct_field( matvar_t *eeg, const char *name, int struct_array_index ){
  matvar_t *tmp; 

  tmp = Mat_VarGetStructField( eeg, (char*)name, BY_NAME, struct_array_index );
  if( tmp->rank != 1 && tmp->dims[0]<1 ){
	 errprintf("field '%s' wrong, rank=%i,tmp->dims[0]=%i\n",name, tmp->rank,tmp->dims[0] );
	 return -1;
  }  
  dprintf( "found: %s=%f\n", name, (((double*)tmp->data)[0]) );
  return (((double*)tmp->data)[0]);
}


/* ------------------------------- WRITER ----------------------------- */



/** \brief write EEG-struct to EEGlab-compatible MATLAB file.

	 This write uses MatIO to create an EEGlab file.
	 It was developed with EEGlab version 6.01b.
	 
	 \todo NOT IMPLEMENTED!

	 \param eeg the struct
	 \param file eeglab .set file
	 \return error code
*/
int write_eeglab_file( EEG* eeg, const char *file ){
  mat_t *mfile;
  matvar_t *meeg=NULL;
  if( !(mfile=Mat_Create( file, NULL )) ){
	 errprintf("Could not open '%s' for writing\n", file);
	 return -1;
  }

  //int dims[2]={1,1}, rank=2;
  //meeg = Mat_VarCreate( "EEG", MAT_T_STRUCT, , rank, dims, NULL, 0 );
  
  Mat_VarWrite( mfile, meeg, 0 );
  return 0;
}

/** \brief save Array-struct as MATLAB .mat file.

	 Currently, the function converts all arrays to
	 double-arrays. Should be easy to implement for arbitrary types if
	 needed.

	 \param a the Array
	 \param varname MATLAB-name of the array
	 \param file the filename
	 \param append if TRUE, the function tries to open() the file, else it is overwritten
	 \return error code
 */
int write_array_matlab( const Array *a, const char *varname, const char *file, bool append ){
  mat_t *mfile;
  matvar_t *marr=NULL;
  int i;
  
  if( append ){
	 if( !(mfile=Mat_Open( file, MAT_ACC_RDWR )) ){
		errprintf("Could not open '%s', creating new file\n", file);
	 } 
  }
  
  if( !mfile ){
	 if( !(mfile=Mat_Create( file, NULL )) ){
		errprintf("Could not open '%s' for writing\n", file);
		return -1;
	 }
  }

  int ndim = MAX(2, a->ndim);
  int *size = (int*)malloc( ndim*sizeof(int));
  if( a->ndim==1 ){
	 size[0]=1;
	 size[1]=a->ndim;
  } else {
	 memcpy( size, a->size, ndim*sizeof(int));
  }
  dprintf("Writing to file '%s' variable '%s', ndim=%i\n", file, varname, ndim);

  /* convert to column major for MATLAB */
  Array *b=array_convert_rowcolmajor( (Array*)a, TRUE );
  
  /* up-cast to DOUBLE - copy of b */
  Array *c=array_new( DOUBLE, b->ndim, b->size );
  for( i=0; i<array_NUMEL(b); i++ ){
	 array_dtype_to_double( array_INDEXMEM1(c,i), array_INDEXMEM1(b,i), b->dtype );
  }
  array_free( b );

  marr = Mat_VarCreate( varname, MAT_C_DOUBLE, MAT_T_DOUBLE, 
								ndim, size, c->data, MEM_CONSERVE /* Array remains owner */
								);
  
  dprintf("mfile=%p, marr=%p\n", mfile, marr );
  int succ=Mat_VarWrite( mfile, marr, 0 );
  dprintf("done writing with succ=%i\n", succ);

  Mat_Close( mfile );
  Mat_VarFree( marr );
  array_free( c );
  free( size );

  return 0;
}
#endif /* MATIO */


