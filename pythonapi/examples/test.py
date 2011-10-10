from pyeegtools import *
from pylab import *

EEG = read_eeglab_file("../install/share/vp19_tt.set");
x = EEG.get_data();
z = x[1,1,];

plot(z, 'b');
running_median(z, 10)

moving_average( z, 20 );
plot( z, 'r');

show();
