/***************************************************************************
 *   Copyright (C) 2010 by Matthias Ihrke                                  *
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
 *   aint with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/**\file io_wav.h
 \brief \ref status_inprogress io_wav WAV file input-output for (EEG)-data.

 The WAV-file format allows multi-channel output.
 Writing the EEG-data to WAV is a convenient way for visualizing
 the data in an audio-visualizer software (e.g. Audacity).
	
 For writing the WAV file, some parameters need to be passed to the
 \ref WavFile struct:

 
 \todo implement matrix-to-wavfile conversion

 */
#ifndef IO_WAV_H
# define IO_WAV_H

#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif
  
 
  
  typedef struct wavfile_t {
	 uint16_t numchannels; /**< Mono = 1, Stereo = 2, etc. */
	 uint32_t samplerate;  
	 uint32_t byterate;    /**<  == SampleRate * NumChannels * BitsPerSample/8 */
	 uint16_t blockalign;  /**< == NumChannels * BitsPerSample/8
									  The number of bytes for one sample including
									  all channels. */
	 uint16_t bits_per_sample;
	 
	 uint32_t data_size;  /**< == NumSamples * NumChannels * BitsPerSample/8
									 This is the number of bytes in the data.*/
	 void *data;          /**< 8-bit samples are stored as unsigned bytes,
									 ranging from 0 to 255. 16-bit samples are stored
									 as 2's-complement signed integers, ranging
									 from -32768 to 32767. */
  } WavFile;
  
  /* -------------- FUNCTIONS ---------------- */
  
  WavFile* wavfile_read ( FILE *f );
  int      wavfile_write( WavFile *w, const char *fname );
  
  void     wavfile_print( FILE *out, WavFile *w );
  void     wavfile_free ( WavFile *w );
  
#ifdef __cplusplus
}
#endif

#endif /* IO_WAV_H */
