[ distmat ] = signaldistmatrix( signal1, signal2, varargin );
%SIGNALDISTMATRIX
%   It calculates the pointwise distance matrix:  d(i,j) = d( s1_i, s2_j )
%
%   [ distmat ] = signaldistmatrix( signal1, signal2 );
%    Calculate a euclidean distance matrix.
%
%   [ distmat ] = signaldistmatrix( signal1, signal2, metric, ... );
%    Calculate a distance matrix with another metric.
%    Valid values for metric are:
%      "euclidean", ...
%    You can pass additional parameters for the distance function along.
%
%Example:
%
%References:
%
%See Also: SIGNALDISTMATRIX
