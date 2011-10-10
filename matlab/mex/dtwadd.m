 [ average ] = dtwadd( s1, s2, varargin );
% DTWADD
%    This function adds two signals either directly computing
%    the warping, or using a previously calculated warping-path.
%
%    [ average ] = dtwadd( s1, s2 );
%        s1 and s2 are matrices or vectors giving the signals to be warped and
%        added.
%
%    [ average ] = dtwadd( s1, s2, path );
%        s1 and s2 are matrices or vectors giving the signals to be
%        added. path is a previously calculated warping path, i.e.
%        a 2 x M integer array giving corresponding points in the
%        signals.
