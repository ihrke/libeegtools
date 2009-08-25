 /* put some more python functionality here */
%pythoncode
%{

import pylab
import math


def plot_eegdata(eeg, plotfilename=None, times=None, channels=None, format="b-"):
    """
    """
    if not channels:
        channels = range(0,eeg.nbchan);
    d = eeg.get_data();
    nc = int(math.ceil(math.sqrt( len(channels) )));
    idx = 1;
    for c in channels:
        pylab.subplot(nc,nc,idx);
        if times:
            pylab.plot(times, d[c], format);
        else:
            pylab.plot(d[c], format);
        pylab.title("Channel %i"%c);
        idx = idx+1;
    if( plotfilename ):
        pylab.savefig( plotfilename );

def plot_eeg(eeg, plotfilename=None, times=None, channels=None, trials=None, format="b-"):
    """
    plotting for EEGdata and EEGdata_trials.
    pylab based.
    """
    if type(eeg).__name__=='EEGdata':
        print "This is an EEGdata struct";
        plot_eegdata( eeg, plotfilename=plotfilename,
                      times=times, channels=channels, format=format );
    elif type(eeg).__name__=="EEGdata_trials":
        print "This is an EEGdata_trials struct"
        
    else:
        print "Plotting does not support '%s'"%(type(eeg).__name__)


        
%}
