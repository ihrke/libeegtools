from scipy import *

double64 = dtype('f8'); # floating point (double64)

def read_raw( filename ):   
    raw = fromfile(file=filename, dtype=double64);
    nchan   = raw[0];
    ntrials = raw[1];
    nsample = raw[2];
    nmarkers= raw[3];
    offset = 4;
    times   = raw[offset:offset+nsample];
    offset += nsample;
    markers = reshape( raw[offset:offset+(ntrials*nmarkers)], (-1, 2) );
    offset += (ntrials*nmarkers);
    data    = reshape( raw[offset:], (nchan, ntrials, nsample) );
    
    return (times,markers,data);

def rmse( realsig, sig):
    """
    root mean-square error
    """
    # r = sqrt( 1/length(realsig) * sum( (realsig-sig).^2));
    return sqrt(1.0/size(sig)*sum(pow(realsig-sig, 2)));

def snr( realsig, sig ):
    """
    Signal-to-noise ratio
    """
    #    s = 10*log10(sum(realsig.^2) / sum((realsig-sig).^2));
    return 10*log10(sum(pow(realsig,2))/(sum( pow(realsig-sig, 2) )));
