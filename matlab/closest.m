function [ index ] = closest( v, x )
%closest return index of closest entry in vector v to scalar x (abs)

[tmp index] = min(abs(v-x));
