% % pop_warpavg() - Timewarped-Average
%
% Usage:
%   >>  OUTEEG = pop_warpavg( INEEG, 'key', 'value',...);
%
% Inputs:
%   INEEG   - input EEG dataset, must be epoched data
%
% Opts:
%   GLOBAL OPTIONS:
%   --------------_
%    chans        - select which channels are to be corrected; []=all
%    rtrigger     - ['trigger1' 'trigger2' ... ] list of triggers that
%                   can be considered as response (if a trial does not
%                   contain response information, an error is raised);
%                   If the list is empty, any trigger following the
%                   0-trigger (stimulus) is considered as a response.
%
%   WAVELET-OPTIONS:
%   ----------------
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
%   TIMEWARP-OPTIONS:
%   -----------------
%    theta1,theta2 - weights for the timewarping norm;
%                    theta1 - weight for amplitude of signal
%                    theta2 - weight for gradient of signal
%    tau           - time in ms to warp before stimulus onset
%    b             - time in ms to warp after response trigger
%    maxiterations - maximum number of iterations (def: 10)
%    ccrit         - convergence-criterion: percentage that the RMSE 
%                    must have changed in order to continue the
%                    iteration (default: 5)
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
% Revision 1.1  2008/03/14 13:57:45  mihrke
% Initial revision
%
%
%

function [EEG, com] = pop_warpavg( EEG, varargin );

% the command output is a hidden output that does not have to
% be described in the header

com = ''; % this initialization ensure that the function will return something
          % if the user press the cancel button            

% display help if not enough arguments
% ------------------------------------
if nargin < 1
	help pop_warpavg;
	return;
end;	

% cmd-line parsing
args = varargin;
% create structure
% ----------------
if ~isempty(args)
   try, g = struct(args{:});
   catch, disp('pop_warpavg(): wrong syntax in function arguments'); return; end;
else
    g = [];
end;

% test the presence of variables
% ------------------------------
try, g.chans;              catch, g.chans=1:EEG.nbchan; end;
try, g.do_wavelet;          catch, g.do_wavelet=1; end;
try, g.L;                  catch, g.L=6; end;
try, g.lambda;             catch, g.lambda='heursure'; end;
try, g.eta;                catch, g.eta='soft'; end;
try, g.sigext;             catch, g.sigext='sym'; end;
try, g.theta1;             catch, g.theta1=1.0; end;
try, g.theta2;             catch, g.theta2=1.0; end;
try, g.rtrigger;           catch, g.rtrigger={}; end;
try, g.tau;                catch, g.tau=100; end;
try, g.b;                  catch, g.b=100; end;
try, g.maxiterations;      catch, g.maxiterations=10; end;
try, g.ccrit;              catch, g.ccrit=5.0; end;

assert((EEG.trials>1 && ~isempty(EEG.epoch)),[' Error: pop_warpavg(): Should be ' ...
                    'applied to epoched data.']);
if nargin <= 1
  % popup window parameters
  % -----------------------
  % -----------------
  uilist = { ...
      { 'style' 'text' 'string' 'General Options', 'fontweight', 'bold'} ...
      { 'style' 'text' 'string' 'Channels ([]=all):' } ...
      { 'style' 'edit' 'string' '' } ...        
      { 'style' 'text' 'string' 'Response Triggers ({}=any in segment):' } ...
      { 'style' 'edit' 'string' '' } ...     
      { 'style' 'text' 'string' 'Wavelet-Options', 'fontweight', 'bold'} ...
		{ 'style' 'text' 'string' 'Perfom Thresholding?'} ...
		{ 'style' 'checkbox' 'string' '' 'value' 1} ...
      { 'style' 'text' 'string' 'L:' } ...
      { 'style' 'edit' 'string' '6' } ...
      { 'style' 'text' 'string' 'lambda:' } ...
      { 'style' 'edit' 'string' 'heursure' } ...
      { 'style' 'text' 'string' 'eta:' } ...
      { 'style' 'edit' 'string' 'soft' } ...
      { 'style' 'text' 'string' 'sigext:' } ...
      { 'style' 'edit' 'string' 'sym' } ...
      { 'style' 'text' 'string' 'Timewarp-Options', 'fontweight', 'bold'} ...
      { 'style' 'text' 'string' 'tau (ms):' } ...
      { 'style' 'edit' 'string' '100.0' } ...
      { 'style' 'text' 'string' 'b (ms):' } ...
      { 'style' 'edit' 'string' '100.0' } ...    
      { 'style' 'text' 'string' 'Max. Iterations:' } ...
      { 'style' 'edit' 'string' '10' } ...
      { 'style' 'text' 'string' 'Convergence Crit. (%):' } ...
      { 'style' 'edit' 'string' '5' } ...
      { 'style' 'text' 'string' 'theta1:' } ...
      { 'style' 'edit' 'string' '1.0' } ...
      { 'style' 'text' 'string' 'theta2:' } ...
      { 'style' 'edit' 'string' '1.0' } ...
           };
  geometry = [1 2 2 1 2 2 2 2 2 1 2 2 2 2 2 2];
  
  result = inputgui( 'geometry', geometry, 'uilist', uilist, 'title', 'Timewarped Averaging  -- pop_warpavg()', ...
                     'helpcom', 'pophelp(''pop_warpavg'')');
  
  chans = eval(sprintf('[ %s ]', result{1}));
  if(~isempty(chans))  g.chans=chans;  end;
  rtrig = eval(sprintf('{ %s }', result{2}));
  if(~isempty(rtrig))  g.rtrigger=rtrig;  end;
  g.do_wavelet = (result{3});
  g.L = eval(result{4});
  g.lambda = (result{5});
  g.eta = (result{6});
  g.sigext = (result{7});
  g.tau =  eval(result{8});
  g.b =  eval(result{9});
  g.maxiterations =  eval(result{10});
  g.ccrit = eval(result{11});
  g.theta1 = eval(result{12});
  g.theta2 = eval(result{13});
end;

l.conventional=0;
l.ti=1;
l.sureshrink=2;
l.heursure=3;
strisin(fieldnames(l), g.lambda)
assert(logical(strisin(fieldnames(l), g.lambda)),[' Error: pop_warpavg(): unknown lambda']);
g.lambda = (l.(g.lambda));

e.hard=0;
e.soft=1;
assert(logical(strisin(fieldnames(e), g.eta)),[' Error: pop_warpavg(): unknown eta']);
g.eta = (e.(g.eta));

s.zeros=0;
s.zerosr=1;
s.sym=2;
s.smooth=3;
assert(logical(strisin(fieldnames(s), g.sigext)),[' Error: pop_warpavg(): unknown sigext']);
g.sigext = (s.(g.sigext));


disp(g)
% Find the RTs from the rtrigger-cell-array
rts = [];
for i=1:EEG.trials
  lat = [EEG.epoch(i).eventlatency{:}];
  trigs = EEG.epoch(i).eventtype;
  index = find(lat==0);
  assert(~isempty(index), sprintf([' Error: pop_warpavg(): no stimulus ' ...
                      'time in epoch %i'],i));
  assert(size(lat,2)>index, sprintf([' Error: pop_warpavg(): no response ' ...
                      'time in epoch %i'],i));
  j = 1;
  if(~isempty(g.rtrigger))
    while(~strisin(g.rtrigger, trigs(index+j)))
      assert(size(lat,2)>=(index+j), sprintf([' Error: pop_warpavg(): no response ' ...
                          'time in epoch %i'],i));
      j=j+1;
    end;
  end;
  rts = [rts lat(index+j)];
end;

%disp(rts);
EEG.warpavg=[];

stop=false;
h = waitbar(0.0,'Please wait...');
for c = g.chans
	if(g.do_wavelet)
		EEG = pop_waveden( EEG, 'chans', c, 'L', g.L, 'sigext', g.sigext, 'lambda', g.lambda, 'eta', g.eta );
	end;
	EEG.warpavg(c,:) = ml_warpavg(EEG.times, reshape(double(EEG.data(c,:,:)), [EEG.pnts ...
                      EEG.trials]), rts, g.L, g.lambda, g.eta, g.sigext, ...
                                g.theta1, g.theta2, g.tau, g.b, g.maxiterations);  
  waitbar(c/(g.chans(end)),h);
end;


if exist('h') delete(h); end

% return the string command
% -------------------------
com='';
%com = sprintf(['EEG = pop_warpavg( EEG, ''chans'', [%s], ''L'', %d, ''lambda'', ''%s'', ''eta'', ''%s'', ''sigext'', ''%s'', ''rtrigger'', ''%s'', ''theta1'', %f, ''theta2'', %f);',...
%               num2str(g.chans), g.L, g.lambda, g.eta, g.sigext], g.rtrigger, ...
%              g.theta1, g.theta2);

return;