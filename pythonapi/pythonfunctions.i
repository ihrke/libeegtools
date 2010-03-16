 /* put more python functionality here */
%pythoncode
%{

import pylab
import math


def plot_eeg(eeg, plotfilename=None, channels=None, trials=None, format="b-"):
    """
    Plotting an EEG struct from pyeegtools.
    pylab based.
    """
    if type(eeg).__name__!='EEG':
        print "This is not an EEG struct";
        return
    
    if not channels:
        channels = range(0,eeg.nbchan);
    d = eeg.get_data();
    times = eeg.get_times();
    nc = int(math.ceil(math.sqrt( len(channels) )));
    idx = 1;
    for c in channels:
        pylab.subplot(nc,nc,idx);
        if times!=None:
            pylab.plot(times, d[c,trials], format);
        else:
            pylab.plot(d[c,trials], format);
        pylab.title("Channel %i"%c);
        idx = idx+1;
    if( plotfilename ):
        pylab.savefig( plotfilename );
    
%}
