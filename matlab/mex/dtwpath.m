 [ path, cumulated_distmat ] = dtwpath( in1, varargin );
%DTWPATH
%   The script can be used to calculate the warping-path for
%   either a given distance matrix between individual points in
%   two signals, or by giving the (potential multidimensional) signals
%   directly.
%
%   It returns the warping path (that is corresponding elements in
%   both series as a 2 x M integer-array).
%
%   [ path, cumulated_distmat] = dtwpath( s1, s2 );
%    s1 and s2 are matrices or vectors giving the signals to be warped.
%
%   [ path, cumulated_distmat ] = dtwpath( s1, s2, metric );
%    s1 and s2 are matrices or vectors giving the signals to be warped.
%    metric is a string giving the distance metric. \todo implement this!
%
%   [ path, cumulated_distmat ] = dtwpath( distmat );
%     distmat is the n x n pointwise distance matrix between signal points.
