from libeegtools import *
import pylab

print "Exemplary mathematical functions"
x=bresenham(0,10, 10,20,intArray(100));

print "Read EEG data from file"
eeg = read_eegtrials_from_raw( "/home/ihrke/sloc/data/test_ttjj.raw");#artdata1.raw" );
print eeg
x = eeg.get_trial(1); # EEGdata class
print x
z = x.get_chan(0); # list containing data from channel 0
d = x.get_data();  # list of lists of all channels


plot_eeg( x, channels=[3,4] );


y = wavelet_denoising( x );
plot_eeg( y, channels=[3,4], format="r-" );
pylab.show();

# pylab.clf()
# y = x.clone();
# eeg_filter_running_median(y, 50);
# plot_eeg( y, channels=[3,4], plotfilename="test2.pdf" );

# pylab.clf()
# y = x.clone();
# plot_eeg( y, channels=[3,4], format="r-");
# eeg_filter_fidlib(y, eeg.sampling_rate, "BpBu4/0.1-80");
# plot_eeg( y, channels=[3,4], plotfilename="test3.pdf" );


# print "Write EEG data to file"


# print "Creating new EEG data"
# x=EEGdata(10,100, 2);

# x=EEGdata_trials(10, 2, 64, 1000);
