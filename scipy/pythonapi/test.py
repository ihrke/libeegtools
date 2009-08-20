from libeegtools import *

x = Complex(1,2);
print x


x=bresenham(0,10, 10,20,intArray(100));

print "mad=%f"%mad( [1,2,3,4,100],5);

x=EEGdata(10,100, 2);
print x;

x=EEGdata_trials(10, 2, 64, 1000);
print x;
