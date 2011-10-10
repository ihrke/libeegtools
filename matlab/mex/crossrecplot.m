[ R ] = crossrecplot( s1, s2, varargin );
%CROSSRECPLOT
%   This script calculates the cross-recurrence plot of two (multi-dimensional)
%   signals.
%
%   The MATLAB-interface is as follows:
%
%   [ R ] = crossrecplot( s1, s2 );
%     s1 and s2 are matrices or vectors giving the signals.
%
%   [ R ] = crossrecplot( s1, s2, epsilon );
%     s1 and s2 are matrices or vectors giving the signals, epsilon is
%     the neighbourhood of the signals. The signals should be normalized before.
%
%   [ R ] = crossrecplot( s1, s2, epsilon, neighbourhood_crit );
%     s1 and s2 are matrices or vectors giving the signals, epsilon is
%     the neighbourhood of the signals or a fixed amount of neighbours (FAN) that
%     is used to calculate the recurrences. Whether an epsilon-ball or a FAN is
%     used depends on neighbourhood_crit which is one of 'fan' or 'epsilon'.
%
% See also: LOSCALC
%
%
% Examples:
% a = sin((1:1000) * 2 * pi/67);
% R=crossrecplot(a, a);
% imagesc( R );
