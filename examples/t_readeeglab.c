/** \example t_readeeglab.c
 *
 * \brief Testing EEGlab-reader
 *
 * Compilation:
 *\code
 *\endcode
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memcpy */
#include "reader.h"
#include "definitions.h"
#include "helper.h"
#include "eeg.h"
#include "wavelet.h"


void plot_eeg( EEG *eeg, const char *format );
/* ---------------------------------------------------------------------------- 
   -- main routine                                                           -- 
   ---------------------------------------------------------------------------- */
int main(int argc, char **argv){
  FILE *f;
  int i, j;
  EEG *eeg, *eeg2, *avg;

  oprintf("Reading file\n");
  //  const char *file="/data/voicekey/eeglab/filtered05_20Hz/eeg081206_1_resampled500hz_filtered.set";
  //  const char *file="/data/eeglabsets_exp6/young/eegtast_01_exp_resampled500hz_filtered.set";
  //  const char *file="/data/voicekey/eeglab/final/eeg061206_5_resampled500hz_filtered_TT.set";
  //  const char *file="/data/exp6test.set";
  //  const char *file="/data/voicekeytest.set";
  const char *file="/media/scratch.local/ihrke/data/vp19_tt.set";

  eeg=read_eeglab_file( file ); 
  eeg_print( stdout, eeg, 3 );
  //  plot_eeg( eeg, "-r");
  //plot_show();


  oprintf("Cloning eeg\n");
  eeg2 = eeg_clone( eeg, EEG_CLONE_ALL );//NODATA|EEG_CLONE_NOMARKERS|EEG_CLONE_NOCHANINFO );
  //  plot_eeg( eeg2, "-r");
  //  plot_show();

  eeg_remove_baseline( eeg2, -500, 0, FALSE );
  //  plot_eeg( eeg2, "-r");
  //  plot_show();

  oprintf("Gibbons Trials\n");
  avg = eeg_gibbons( eeg, 0, 1, 1.0 );
  eeg_print( stdout, avg, 3 );
  plot_eeg( avg, "-r");

#ifdef PLOTTER
  plot_show();
#endif

#if 0
  oprintf("Extract Trials\n");
  int trials[10] = {1,2};
  eeg2 = eeg_extract_trials( eeg2, trials, 2, FALSE );
  eeg_print( stdout, eeg2, 3 );


  oprintf("Extract channels\n");
  int chans[10]={1,2,5,65};
  eeg2 = eeg_extract_channels( eeg2, chans, 1, FALSE );
  plot_eeg( eeg2, "-r");
  plot_show();

  oprintf("Filtering\n");
  //  eeg2 = eeg_filter_fidlib( eeg2, "BpBu4/0.5-20", FALSE );
  WaveletParameters P;
  P = wavelet_init();
  P.threshselfct = conventional_thresholding;
  eeg2 = eeg_wavelet_denoise( eeg2, P, FALSE );
  eeg_print( stdout, eeg2, 3 );
  plot_eeg( eeg2, "-r");
  plot_show();
#endif

  eeg_free( eeg );
  eeg_free( eeg2 );


  oprintf("Read eeg from raw\n");
  //  const char *rawfile="/media/scratch.local/ihrke/data/artdata1.raw";
  const char *rawfile="/data/artdat.raw";
  eeg = read_eeg_from_raw( rawfile );
  eeg_print( stdout, eeg, 3 );
  return 0;
}
void plot_eeg( EEG *eeg, const char *format ){
#ifdef PLOTTER
  int chan, c, sb;
  double sbpos[4], sbsize[4];
  double psize = 10;
  double position[2], size[2]={10,10};
  int nm, i;  
  double xlinex[2], yliney[2], xliney[2]={0};
	 
  xlinex[0] = (eeg->times)?eeg->times[0]:0;
  xlinex[1] = (eeg->times)?eeg->times[eeg->n-1]:eeg->n-1;
  yliney[0] = -100; yliney[1] = 100;


  for( c=0; c<eeg->nbchan; c++ ){
	 chan = c;
	 position[0] = eeg->chaninfo[chan].y; 
	 position[1] = eeg->chaninfo[chan].x; 

	 sbpos[0] = position[0]; sbpos[1]=position[1];
	 sbpos[2] = sbpos[0]+psize;  sbpos[3]=sbpos[1]+psize;
	 sbsize[0]= (eeg->times)?eeg->times[0]:0; 
	 sbsize[1]= (eeg->times)?eeg->times[eeg->n-1]:eeg->n-1;//eeg->n;
	 sbsize[2]=-100; sbsize[3]=100;


	 sb = plot_subplot_create(sbpos, sbsize);
	 plot_subplot_select(sb); 
	 plot_format( xlinex, xliney, 2, "k-" );
	 plot_format( xliney, yliney, 2, "k-" );
	 plot_format_nocopy( eeg->times, eeg->data[chan][0], eeg->n, format );

	 /* nm = eeg->nmarkers; */
	 /* for( i=0; i<nm; i++ ){ */
	 /* 	/\* dprintf(" i=%i, (x,y)=(%*\/ */
	 /* 	plot_format_nocopy( &(eeg->times[eeg->markers[i]]), &(eeg->d[chan][eeg->markers[i]]), 1, "kO" ); */
	 /* } */
	 plot_subplot_select(-1);
  }
#endif
}
