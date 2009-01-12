% % pop_waveden() - Wavelet-based Denoising of Epoched signals.
%
% Usage:
%   >>  OUTEEG = pop_waveden( INEEG, 'key', 'value',...);
%
% Inputs:
%   INEEG   - input EEG dataset, must be epoched data
%
% Opts:
%    chans        - select which channels are to be corrected; []=all
%    L            - first wavelet-level to threshold; default=6;
%    lambda       - threshold-selection scheme;
%                   [conventional|ti|sureshrink|heursure(**)]
%                   conventional: lambda = sigma*sqrt(2\log_e{N})
%                   ti: lambda = \sigma \sqrt{2\log_e{(N\log_2{N})}}
%                   sureshrink: adaptive estimation of optimal threshold
%                        using SURE 
%                   heursure: heuristic implementation of sureshrink
%                        with fallback-mechanism to conventional
%   eta           - thresholding function; [hard|soft(**)]
%   sigext        - signal-extension scheme; \
%                   [zeros|zerosr|sym(**)|smooth]
%                   zeros:  [1 2 3 - - -] -> [0 1 2 3 0 0] 
%                   zerosr: [1 2 3 - - -] -> [1 2 3 0 0 0]
%                   sym: [1 2 3 - - - - -] -> [2 1 1 2 3 3 2 1]
%                   smooth: [1 2 3 - - - - -] -> [1 1 1 2 3 3 3 3]
%
%
% Outputs:
%   OUTEEG  - output dataset
%
% See also:
%   EEGLAB 

% Copyright (C) 2007  Matthias Ihrke <mihrke@uni-goettingen.de>
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% $Log$
% Revision 1.1  2009/01/12 11:34:20  mihrke
% documentation,
% rearrangements,
% eeglab
%
% Revision 1.1.1.1  2008/03/14 13:57:45  mihrke
% Initial Import
%
%
%

function [EEG, com] = pop_waveden( EEG, varargin );

% the command output is a hidden output that does not have to
% be described in the header

com = ''; % this initialization ensure that the function will return something
          % if the user press the cancel button            

% display help if not enough arguments
% ------------------------------------
if nargin < 1
	help pop_waveden;
	return;
end;	

% cmd-line parsing
args = varargin;
% create structure
% ----------------
if ~isempty(args)
   try, g = struct(args{:});
   catch, disp('pop_waveden(): wrong syntax in function arguments'); return; end;
else
    g = [];
end;

% test the presence of variables
% ------------------------------
try, g.chans;              catch, g.chans=1:EEG.nbchan; end;
try, g.L;                  catch, g.L=6; end;
try, g.lambda;             catch, g.lambda='heursure'; end;
try, g.eta;                catch, g.eta='soft'; end;
try, g.sigext;             catch, g.sigext='sym'; end;

assert((EEG.trials>1),[' Error: pop_waveden(): Denoising should be ' ...
                    'applied to epoched data.']);

if nargin <= 1
  % popup window parameters
  % -----------------------
  % -----------------
  uilist = { ...

      { 'style' 'text' 'string' 'Channels ([]=all):' } ...
      { 'style' 'edit' 'string' '' } ...        
      { 'style' 'text' 'string' 'L:' } ...
      { 'style' 'edit' 'string' '6' } ...
      { 'style' 'text' 'string' 'lambda:' } ...
      { 'style' 'edit' 'string' 'heursure' } ...
      { 'style' 'text' 'string' 'eta:' } ...
      { 'style' 'edit' 'string' 'soft' } ...
      { 'style' 'text' 'string' 'sigext:' } ...
      { 'style' 'edit' 'string' 'sym' } ...
           };
  geometry = [ 2 2 2 2 2 ];
  
  result = inputgui( 'geometry', geometry, 'uilist', uilist, 'title', 'Wavelet-based Denoising -- pop_waveden()', ...
                     'helpcom', 'pophelp(''pop_waveden'')');
  
  chans = eval(sprintf('[ %s ]', result{1}));
  if(~isempty(chans))
    g.chans=chans;
  end;
  g.L = eval(result{2});
  g.lambda = (result{3});
  g.eta = (result{4});
  g.sigext = (result{5});
end;

l.conventional=0;
l.ti=1;
l.sureshrink=2;
l.heursure=3;
strisin(fieldnames(l), g.lambda)
assert(logical(strisin(fieldnames(l), g.lambda)),[' Error: pop_waveden(): unknown lambda']);
g.lambda = (l.(g.lambda));

e.hard=0;
e.soft=1;
assert(logical(strisin(fieldnames(e), g.eta)),[' Error: pop_waveden(): unknown eta']);
g.eta = (e.(g.eta));

s.zeros=0;
s.zerosr=1;
s.sym=2;
s.smooth=3;
assert(logical(strisin(fieldnames(s), g.sigext)),[' Error: pop_waveden(): unknown sigext']);
g.sigext = (s.(g.sigext));



stop=false;
h = waitbar(0.0,'Please wait...');
for c = g.chans
  EEG.data(c,:,:) = reshape(ml_denoise(reshape(double(EEG.data(c,:,:)), [EEG.pnts ...
                   EEG.trials]), g.L, g.lambda, g.eta, g.sigext), [1 ...
                      EEG.pnts EEG.trials]);
  
  waitbar(c/(g.chans(end)),h);
end;


if exist('h') delete(h); end

% return the string command
% -------------------------
com = sprintf(['EEG = pop_waveden( EEG, ''chans'', [%s], ''L'', %d, ''lambda'', ''%s'', ''eta'', ''%s'', ''sigext'', ''%s'');',...
               num2str(g.chans), g.L, g.lambda, g.eta, g.sigext]);

return;