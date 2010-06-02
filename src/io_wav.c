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

#include "io_wav.h"
#include <stdlib.h>

/** \brief list content of WavFile struct.
	 \param out the output stream (e.g. stderr)
	 \param w the wavfile struct
 */
void wavfile_print( FILE *out, WavFile *w ){
  fprintf( out, "WAV '%p'\n"
			  " numchannels     = %d\n"
			  " samplerate      = %d\n"
			  " byterate        = %d\n"
			  " blockalign      = %d\n"
			  " bits_per_sample = %d\n"
			  " data_size       = %d\n"
			  " data            = %p\n",
			  w, w->numchannels, w->samplerate, w->byterate, 
			  w->blockalign, w->bits_per_sample, w->data_size,
			  w->data );
}

/** \brief free WavFile struct and data.
	 \param w the WavFile struct
 */
void     wavfile_free ( WavFile *w ){
  if( w ){
	 if( w->data )
		free( w->data );
	 free( w );
  }
}


/* ------------- READER --------------------- */

/** \brief read WAV-File.

	see https://ccrma.stanford.edu/courses/422/projects/WaveFormat/
	for details on the wav-file format
	\param f a FILE pointer
	\return the wavfile (freshly allocated) or NULL (failure)
 */
WavFile* wavfile_read( FILE *f ){
  uint32_t buf;
  WavFile *w;
  int stopflag=0;

  if( ffread( &buf, 1, 4, f ) &&  buf!=0x46464952 ){ /* "RIFF" */
	 stopflag=1;
  }
  ffread( &buf, 1, 4, f ); /* trash */
  if( ffread( &buf, 1, 4, f ) &&  buf!=0x45564157 ){ /* 'WAVE' */
	 stopflag=1;
  }
  if( ffread( &buf, 1, 4, f ) &&  buf!=0x20746d66 ){ /* 'fmt ' */
	 stopflag=1;
  }
  
  if( stopflag ){
	 fprintf( stderr, "No WAV-File\n" );
	 return NULL;
  }
  if( ffread( &buf, 1, 4, f ) && buf!=16 ){ /* no compression */
	 fprintf( stderr, "Cannot handle compression...\n" );
	 return NULL;
  }

  w = (WavFile*) malloc( sizeof(WavFile) );
  ffread( w, 1, 16, f );

  if( ffread( &buf, 1, 4, f ) && buf!=0x61746164 ){ /* 'data' */
	 fprintf( stderr, "No data\n" );
	 return NULL;
  }

  if( !ffread( &(w->data_size), 1, 4, f ) )
	 return NULL;
  
  w->data = (void*) malloc( w->data_size );
  ffread( w->data, 1, w->data_size, f );

  return w;
}



/* ------------- WRITER --------------------- */

/**\brief write WAV-File.

	see https://ccrma.stanford.edu/courses/422/projects/WaveFormat/
	for details on the wav-file format
	\param w the wavfile
	\param fname the name of the output file
	\return 0 success, else failure
 */
int      wavfile_write( WavFile *w, const char *fname ){
  FILE *f;
  uint32_t buf[5];
  uint16_t buf2;
  if( !(f = fopen( fname, "w" )) ){
	 fprintf( stderr, "Could not open file '%s'\n", fname );
	 return -1;
  }

  buf[0] = 0x46464952; /* RIFF */
  buf[1] = (uint32_t)(36+w->data_size);
  buf[2] = 0x45564157; /* WAVE */
  buf[3] = 0x20746d66; /* fmt */
  buf[4] = 16;         /* (PCM) */
  fwrite( buf, 1, 5*4, f );
  buf2   = 1;
  fwrite( &buf2, 1,   2, f );
  fwrite( &(w->numchannels), 1,  2, f ); 
  fwrite( &(w->samplerate) , 1,  4, f ); 
  fwrite( &(w->byterate),    1,  4, f ); 
  fwrite( &(w->blockalign),  1,  2, f ); 
  fwrite( &(w->bits_per_sample), 1,  2, f ); 
  buf[0] = 0x61746164; /* data */
  fwrite( buf, 1,   4, f );
  fwrite( &(w->data_size), 1, 4, f );
  fwrite( w->data, 1, w->data_size, f );
  fclose( f );

  return 0;
}
