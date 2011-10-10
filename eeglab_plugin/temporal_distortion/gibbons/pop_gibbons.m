% % pop_gibbons() - Calculates Response-Time-Corrected Average as 
%                   detailed in Gibbons+Stahl (2007) 'Response-time
%                   Corrected averaging of event-related potentials'
%                   in Clinical Neurophysiology, 118, 197-208
%
% Usage:
%   >>  OUTEEG = pop_gibbons( INEEG, 'key', 'value',...);
%
% Inputs:
%   INEEG   - input EEG dataset
%
% Opts:
%   response_events              - either a scalar, or a cell-array
%                                  of strings giving the markers that
%                                  indicate that a response happened
%   response_from_stimulus_onset - [0|1]; if set, response_events is
%                                  assumed to be a single number n that
%                                  gives the n'th timemarker after
%                                  stimulus-onset
%   chans                        - array of channels to use
%   k                            - parameter of the method, see the paper.
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
% $Log: pop_gibbons.m,v $
% Revision 1.1  2009/09/02 09:47:22  mihrke
% There is going to be major restructuring of the API.
% Most EEG-functions will be called differently or get different arguments,
% because the main structs are going to change.
% We will move to 0.5 after this.
%
% Revision 1.1  2009/01/12 15:03:53  mihrke
% gibbons eegplugin
%
%

function [EEG, com] = pop_gibbons( EEG, varargin );

% the command output is a hidden output that does not have to
% be described in the header

com = ''; % this initialization ensure that the function will return something
          % if the user press the cancel button            

% display help if not enough arguments
% ------------------------------------
if nargin < 1
	help pop_gibbons;
	return;
end;	

% cmd-line parsing
args = varargin;
% create structure
% ----------------
if ~isempty(args)
   try, g = struct(args{:});
   catch, disp('pop_gibbons(): wrong syntax in function arguments'); return; end;
else
    g = [];
end;

% test the presence of variables
% ------------------------------
try 
  g.response_from_stimulus_onset; 
catch, 
  g.response_from_stimulus_onset=1;
end;
try
  g.response_events;    
catch 
  if g.response_from_stimulus_onset
    g.response_events=1;
  else 
    g.response_events={};
  end;
end;
try, g.chans;              catch, g.chans=1:EEG.nbchan; end;
try, g.k;                  catch, g.k=1; end;
  


assert(~isempty(EEG.epoch),[' Error: pop_gibbons(): Gibbons Response-time ' ...
                    'corrected average can only be applied to epoched data']);

if nargin < 2 
  cbevent = ['if ~isfield(EEG.event, ''type'')' ...
             '   errordlg2(''No type field'');' ...
             'else' ...
             '   if isnumeric(EEG.event(1).type),' ...
             '        [tmps,tmpstr] = pop_chansel(unique([ EEG.event.type ]));' ...
             '   else,' ...
             '        [tmps,tmpstr] = pop_chansel(unique({ EEG.event.type }));' ...
             '   end;' ...
             '   if ~isempty(tmps)' ...
             '       set(findobj(''parent'', gcbf, ''tag'', ''response_events''), ''string'', tmpstr);' ...
             '   end;' ...
             'end;' ...
             'clear tmps tmpv tmpstr tmpfieldnames;' ];

	% popup window parameters
	% -----------------------
    	% -----------------
   uilist = { ...
       { 'style' 'text' 'string' 'Response-Marker(s):' } ...
       { 'style' 'edit' 'string' '' 'tag' 'response_events' } ...  
       { 'style' 'pushbutton' 'string' '...' 'callback' cbevent } ...
       { 'style' 'checkbox' 'string' '' 'value' 1 } ...
       { 'style' 'text' 'string' 'if checked, use event from stimulus onset' } ...
       { 'style' 'text' 'string' '(Response-Marker field is a single number, then)' } ...
       { 'style' 'text' 'string' 'Channels ([]=all):' } ...
       { 'style' 'edit' 'string' '' } ...        
       { 'style' 'text' 'string' 'k:' } ...
       { 'style' 'edit' 'string' '1' } ...
     };
   geometry = { [2 1 0.5] [0.1 1] [1] [2 1] [2 1]};
    
   
   result = inputgui( 'geometry', geometry, 'uilist', uilist, 'title', ...
                      'Response-Time Corrected Averaging (Gibbons et al) -- pop_gibbons()', ...
                      'helpcom', 'pophelp(''pop_gibbons'')');    
   if length(result) == 0, return; end;
   g.response_from_stimulus_onset=(result{2});
   if g.response_from_stimulus_onset
     g.response_events=eval(result{1});
   else
     g.response_events=result{1};
     if isstr(g.response_events)
       x = g.response_events;
       g.response_events={};
       g.response_events{1} = x;
     end;
   end;
   g.chans = eval(sprintf('[ %s ]', result{3}));
   if isempty(g.chans), g.chans = 1:EEG.nbchan; end;
   g.k = eval(result{4});
end;

[electrodes samples trials] = size(EEG.data);
zero_idx = closest( EEG.times, 0 );
trialRT = []; % vector of response times
trialRT_idx = []; % vector of indices for responses for each trial
if g.response_from_stimulus_onset
  for i=1:trials
    zero_marker = find([EEG.epoch(i).eventlatency{:}]==0);    
    resp_marker = zero_marker+g.response_events;
    tmp = EEG.epoch(i).eventlatency(resp_marker);
    trialRT     = [trialRT tmp{:}];
    trialRT_idx = [trialRT_idx closest(EEG.times, tmp{:})];
  end;
else 
  for i=1:trials
    ev = EEG.epoch(i).eventtype(:);
    for s=g.response_events' % find the first match
      resp_marker = find_str_in_cellstr( ev, s );
      if resp_marker~=0
        break;
      end;
    end;
    assert(resp_marker~=0, sprintf(['Error: pop_gibbons(): Response not in ' ...
                        'time-window for trial-No. %i'], i));
    
    tmp = EEG.epoch(i).eventlatency(resp_marker);
    trialRT     = [trialRT tmp{:}];
    trialRT_idx = [trialRT_idx closest(EEG.times, tmp{:})];
  end;
end;
meanRT = mean(trialRT);
meanRT_idx = closest( EEG.times, meanRT );

EEG.data(g.chans,:,:) = gibbons( EEG.data(g.chans,:,:),...
                                 EEG.times, trialRT, g.k );


% return the string command
% -------------------------

return;