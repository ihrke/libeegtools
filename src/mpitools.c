/***************************************************************************
 *   Copyright (C) 2010 by Matthias Ihrke   *
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

#include "mpitools.h"
#include "eeg.h"

ulonglong_t sizeof_channelinfos( const ChannelInfo *chaninfo, int nbchan ){
  if( !chaninfo )
	 return (ulonglong_t)0;
  else 
	 return (ulonglong_t)nbchan*sizeof(ChannelInfo);
}

/** return the number of bytes of the bytestream created by
	 calling eeg_to_stream.
 */
ulonglong_t sizeof_eegstream( const EEG *eeg ){
  ulonglong_t n, tmp;
  int i;

  n=0;
  n += sizeof(unsigned int) + strlen( eeg->filename ); /* length+filename */
  n += sizeof(unsigned int) + strlen( eeg->comment );  /* length+comment */
  n += 3*sizeof(unsigned int); /* nbchan,ntrials,n */
  n += sizeof(double); /* sampling_rate */

  n += 1; /* if 0, no times array is present, else it is */
  if( eeg->times ){
	 n += eeg->n*sizeof(double); /* times-array */
  }
  n += sizeof(ulonglong_t); /* nchannelinfo */  
  tmp = sizeof_channelinfos( eeg->chaninfo, eeg->nbchan );
  if( tmp>0 )
	 n+=tmp;
  /* data */
  n += eeg->ntrials*eeg->nbchan*eeg->n*sizeof(double);
  /* markers */
  n += 1; /* markers present? */
  if( eeg->markers ){
	 n += eeg->ntrials*sizeof(unsigned int); /* nmarkers */
	 for( i=0; i<eeg->ntrials; i++ ){
		n += eeg->nmarkers[i]*sizeof(unsigned int);
		if( eeg->marker_labels ){
		  n += eeg->nmarkers[i]*MAX_LABEL_LENGTH;
		}
	 }
  } 

  return n;
}

/** convert an EEG-struct to a byte-stream. The stream can be used
	 to pass the struct using MPI_Send, e.g.

	 \code
	 char *eegstream;
	 long nstream;
	 eegstream = eeg_to_stream( eeg, NULL, &nstream );
  	 MPI_Send( (void*)(eegstream), nstream, MPI_CHAR, proc,
 	            EEGSEND_TAG, MPI_COMM_WORLD );
	 \endcode

	 \param eeg
	 \param stream either enough memory to hold the stream, or NULL (allocated)
	 \param n is filled with the length of stream
	 \return stream
 */
char* eeg_to_stream( const EEG *eeg, char *stream, ulonglong_t *n ){
  char *s, /* stream ptr */
	 *t;    /* temporary ptr */
  unsigned int t_uint; 
  ulonglong_t t_ulong;
  int c,trial,i;

  s = stream;
  /*-------- calculate size of stream -----------------*/
  *n=sizeof_eegstream( eeg );
  dprintf("stream is n=%llu long (ulonlong)\n", *n );
  
  if( !s ){
	 s = (char*) malloc( *n );
  }

  /*-------fill stream ------------*/

  /* filename */
  t=s;
  t_uint = strlen( eeg->filename );
  memcpy( t, &t_uint, sizeof( unsigned int ) );
  t += sizeof( unsigned int );
  memcpy( t, eeg->filename, t_uint );
  t += t_uint;
  /* comment */ 
  t_uint = strlen( eeg->comment );
  memcpy( t, &t_uint, sizeof( unsigned int ) );
  t += sizeof( unsigned int );
  memcpy( t, eeg->comment, t_uint );
  t += t_uint;
  /* nbchan,ntrials,n */
  memcpy( t, &(eeg->nbchan), sizeof( unsigned int ) );
  t += sizeof(unsigned int);
  memcpy( t, &(eeg->ntrials), sizeof( unsigned int ) );
  t += sizeof(unsigned int);
  memcpy( t, &(eeg->n), sizeof( unsigned int ) );
  t += sizeof(unsigned int);
  /* sampling rate */
  memcpy( t, &(eeg->sampling_rate), sizeof( double ) );
  t += sizeof( double );
  /* times array */
  if( eeg->times ){
	 *t=1; t++;
	 memcpy( t, eeg->times, eeg->n*sizeof(double) );
	 t += eeg->n*sizeof( double );
  } else {
	 *t=0; t++;
  }
  /* nchannelinfo */
  t_ulong = sizeof_channelinfos( eeg->chaninfo, eeg->nbchan );
  memcpy( t, &t_ulong, sizeof( ulonglong_t ) );
  t += sizeof( ulonglong_t );
  /* channelinfos */
  if( t_ulong>0 ){
	 for( i=0; i<eeg->nbchan; i++ ){
		memcpy( t, &(eeg->chaninfo[i]), sizeof( ChannelInfo ) );
		t += sizeof( ChannelInfo );
	 }
  }
  /* data */
  for( c=0; c<eeg->nbchan; c++ ){
	 for( trial=0; trial<eeg->ntrials; trial++ ){
		memcpy( t, eeg->data[c][trial], eeg->n*sizeof(double) );
		t += eeg->n*sizeof( double );
	 }
  }
  /* markers */
  if( eeg->nmarkers ){
	 *t=1; t++;
	 memcpy( t, eeg->nmarkers, eeg->ntrials*sizeof(unsigned int) );
	 t += eeg->ntrials*sizeof(unsigned int);
	 for( i=0; i<eeg->ntrials; i++ ){
		memcpy( t, eeg->markers[i], eeg->nmarkers[i]*sizeof(unsigned int) );
		t += eeg->nmarkers[i]*sizeof(unsigned int);
	 }
	 if( eeg->marker_labels ){
		*t=1; t++;
		for( trial=0; trial<eeg->ntrials; trial++ ){
		  for( i=0; i<eeg->nmarkers[trial]; i++ ){
			 memcpy( t, eeg->marker_labels[trial][i], MAX_LABEL_LENGTH );
			 t += MAX_LABEL_LENGTH;
		  }
		}
	 } else {
		*t=0; t++;
	 }
  } else {
	 *t=0; t++;
  }

  return s;
}
/** convert stream created by eeg_to_stream() back to 
	 an EEG-struct.

	 \code

	 \endcode

	 \param stream 
	 \param n length of stream
	 \param eeg either an appropriate EEG-set, or NULL (allocated)
	 \return eeg-set
 */
EEG*  stream_to_eeg( const char *stream, ulonglong_t n, EEG *eeg ){
  unsigned int nbchan, ntrials, nsamples;
  unsigned int t_uint;
  ulonglong_t t_ulong;
  char *filename;
  char *comment;
  int c,i,j;

  const char *t; /* temporary pointer to stream */

  t=stream;
  /* filename */
  memcpy( &t_uint, t, sizeof(unsigned int) );
  t+=sizeof(unsigned int);
  filename = (char*) malloc( t_uint+1 );
  memcpy( filename, t, t_uint );
  filename[t_uint]='\0';
  t+=t_uint; 
  /* comment */
  memcpy( &t_uint, t, sizeof(unsigned int) );
  t+=sizeof(unsigned int);
  comment = (char*) malloc( t_uint+1 );
  memcpy( comment, t, t_uint );
  comment[t_uint]='\0';
  t+=t_uint; 
  /* nbchan, ntrials,n */
  memcpy( &nbchan, t, sizeof(unsigned int)); 
  t+=sizeof(unsigned int);
  memcpy( &ntrials, t, sizeof(unsigned int)); 
  t+=sizeof(unsigned int);
  memcpy( &nsamples, t, sizeof(unsigned int)); 
  t+=sizeof(unsigned int);
  
  if( !eeg ){
	 eeg = eeg_init( nbchan, ntrials, nsamples );
  } else {
	 eeg->ntrials = ntrials;
	 eeg->nbchan = nbchan;
	 eeg->n = nsamples;
  }
  eeg->filename = filename;
  eeg->comment = comment;
  /* sampling_rate */
  memcpy( &(eeg->sampling_rate), t, sizeof( double ) );
  t += sizeof(double);
  
  /* times array */
  if(*t>0){ /* there is the times array */
	 t++;
	 eeg->times = (double*)malloc( nsamples*sizeof(double) );
	 memcpy( eeg->times, t, eeg->n*sizeof(double) );
	 t+=eeg->n*sizeof(double);
  } else {
	 t++;
	 eeg->times = NULL;
  }
  /* chaninfo */
  memcpy( &t_ulong, t, sizeof( ulonglong_t ) );
  t += sizeof( ulonglong_t );
  if( t_ulong>0 ){
	 eeg->chaninfo = (ChannelInfo*) malloc( eeg->nbchan*sizeof(ChannelInfo) );
	 for( i=0; i<eeg->nbchan; i++ ){
		memcpy( &(eeg->chaninfo[i]), t, sizeof(ChannelInfo) );
		t+=sizeof(ChannelInfo);
	 }
  }
  /* data */
  for( c=0; c<nbchan; c++ ){
	 for( i=0; i<ntrials; i++ ){
		memcpy( eeg->data[c][i], t, nsamples*sizeof(double) );
		t+=nsamples*sizeof(double);
	 }
  }
  /* markers */
  if( *t>0 ){
	 t++;
	 eeg->nmarkers = (unsigned int*)malloc( ntrials*sizeof(unsigned int) );
	 eeg->markers = (unsigned int**)malloc( ntrials*sizeof(unsigned int*));
	 memcpy( eeg->nmarkers, t, ntrials*sizeof(unsigned int) );
	 t += ntrials*sizeof(unsigned int);
	 for( i=0; i<ntrials; i++ ){
		eeg->markers[i] = (unsigned int*) malloc( eeg->nmarkers[i]*sizeof(unsigned int) );
		memcpy( eeg->markers[i], t, eeg->nmarkers[i]*sizeof(unsigned int) );
		t += eeg->nmarkers[i]*sizeof(unsigned int);
	 }
	 /* marker labels */
	 if( *t>0 ){
		t++;
		eeg->marker_labels = (char***)malloc( ntrials*sizeof(char**) );
		for( i=0; i<ntrials; i++ ){
		  eeg->marker_labels[i] = (char**)malloc( eeg->nmarkers[i]*sizeof(char*) );
		  for( j=0; j<eeg->nmarkers[i]; j++ ){
			 eeg->marker_labels[i][j] = (char*) malloc( MAX_LABEL_LENGTH );
			 memcpy( eeg->marker_labels[i][j], t, MAX_LABEL_LENGTH );
			 t+=MAX_LABEL_LENGTH;
		  }
		}
	 }	
  }

  return eeg;
}
