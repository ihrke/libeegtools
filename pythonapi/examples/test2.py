from pylab import *
from pyeegtools import *

x = arange(0,10);# array([1,2,3])
y = arange(10,20);#array([2,3,4])

#x=rand(1000)
print vectordist_euclidean( x, y );
print (x-y)^2
