[ distmat ] = distmatrix( signal, varargin );
%DISTMATRIX
%   It calculates the NxN distance matrix of each point of the
%   (multi-)dimensional signal. The signal is N x p, where N is
%   the number of observations and p the number of dimensions.
%
%   [ distmat ] = distmatrix( signal );
%    Calculate a euclidean distance matrix.
%
%   [ distmat ] = distmatrix( signal, metric, ... );
%    Calculate a distance matrix with another metric.
%    Valid values for metric are:
%      "euclidean", "euclidean_normalized", "dtw" ...
%    You can pass additional parameters for the distance function.
%
%Example:
%
%References:
%
%See Also: SIGNALDISTMATRIX

