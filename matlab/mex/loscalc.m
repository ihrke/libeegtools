path=loscalc( recplot, varargin );
% LOSCALC
%    Calculate the Line-of-synchrony of a cross-recurrence plot.
%
%    los=loscalc( recplot );
%    use the DTW-based approach by Ihrke et al.
%
%    los=loscalc( recplot, 'ihrke2009' );
%    use the DTW-based approach by Ihrke et al.
%
%    los=loscalc( recplot, 'marwan', [dx], [dy] );
%    use the two-step algorithm from Marwan et al.
%    Parameters are optional, default=4
%
% Examples:
% a = sin((1:1000) * 2 * pi/67);
% R=crossrecplot(a, a, 200, 'fan');
% los=loscalc(R);
% imagesc(R);
% axis image; hold on;
% plot( los(1,:), los(2,:), 'w');
%
% x=-10:0.01:10;
% phi=1; psi=1; a=2;
% f=sin(phi*x);
% g=sin(phi*x+a*sin(psi*x));
% R = crossrecplot( f, g, 100, 'fan');
% los=loscalc(R);
% imagesc(R);
% axis image; hold on;
% plot( los(1,:), los(2,:), 'w');
%
%
% See also: CROSSRECPLOT
%
% References:
%
%	Marwan et al. Cross recurrence plot based synchronization
%	of time series. Nonlinear Processes in Geophysics (2002) vol. 9 (3-4) pp. 325-331.
%
%	Matthias Ihrke, Hecke Schrobsdorff and J. Michael Herrmann: Recurrence-Based
%	Synchronization of Single Trials for EEG-Data Analysis. Lecture Notes on Computer
%	Science 5788, Intelligent Data Engineering and Automated Learning - IDEAL 2009.
%	118-125 doi:10.1007/978-3-642-04394-9
