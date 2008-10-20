#!/usr/bin/env python
"""
Extracting information from Raw-Files.
"""
import sys, os, ctypes, struct


raw = sys.argv[1];
sdouble = ctypes.sizeof(ctypes.c_double);

f = open( raw, 'rb' );

numchan=struct.unpack( 'd', f.read( 1*sdouble ))[0];
numtrials=struct.unpack( 'd',  f.read( 1*sdouble ))[0];
numsamples=struct.unpack( 'd', f.read( 1*sdouble ))[0];
nmarkers=struct.unpack( 'd', f.read( 1*sdouble ))[0];
range_start = struct.unpack( 'd', f.read( 1*sdouble ))[0];
step = abs(range_start-struct.unpack( 'd', f.read( 1*sdouble ))[0]);
f.read( int((numsamples-3)*sdouble) );
range_end = struct.unpack( 'd', f.read( 1*sdouble ))[0];
markers = struct.unpack( 'ddddd', f.read( 5*sdouble ));
f.read( int( (nmarkers*numtrials-5)*sdouble ) );
eeg = struct.unpack( 'ddddd', f.read( 5*sdouble ));
f.close();

print "Contents of %s: "%(raw)
print "---------------------------"
print "Number of channels            : %i"%numchan
print "Number of Trials              : %i"%numtrials
print "Number of Samples (per trial) : %i"%numsamples
print "Number of markers (per trial) : %i"%nmarkers
print "Times Array                   : %.2f ... %.2f (step=%.2f)"%(range_start, range_end, step)
print "Sampling Rate                 : %.2f Hz"%(1.0/float(step)*1000)
print "First 5 entries in markers    : %i, %i, %i, %i,%i ..."%(markers[0], markers[1], markers[2],\
                                                        markers[3], markers[4])
print "First 5 entries in EEG-data   : %.2f, %.2f, %.2f, %.2f, %.2f ..."%(eeg[0], eeg[1],\
                                                                          eeg[2], eeg[3], eeg[4])
