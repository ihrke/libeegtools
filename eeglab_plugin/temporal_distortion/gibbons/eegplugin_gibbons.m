% Copyright (C) 2008  Matthias Ihrke <mihrke@uni-goettingen.de>
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
% Revision 1.1  2009/09/02 09:47:22  mihrke
% There is going to be major restructuring of the API.
% Most EEG-functions will be called differently or get different arguments,
% because the main structs are going to change.
% We will move to 0.5 after this.
%
% Revision 1.4  2009/01/12 15:03:53  mihrke
% gibbons eegplugin
%
% Revision 1.3  2009/01/12 11:50:07  mihrke
% adsf
%
% Revision 1.2  2009/01/12 11:43:20  mihrke
% reomve
%
% Revision 1.1  2009/01/12 11:34:20  mihrke
% documentation,
% rearrangements,
% eeglab
%
%
%
function eegplugin_gibbons( fig, try_strings, catch_strings, varargin);


% build command for menu callback
gibbonscmd = [ 'EEG = pop_gibbons(EEG); [ALLEEG EEG CURRENTSET]=eeg_store(ALLEEG,EEG);'...
               'eeglab redraw;'];
fgibbonscmd =  [ try_strings.no_check gibbonscmd ];
fgibbonscmd = [fgibbonscmd  'LASTCOM = ''' gibbonscmd ''';' ];
fgibbonscmd = [fgibbonscmd  catch_strings.store_and_hist ];


% if called from super-plugin, put it under temporal distortion
args = varargin;
if ~isempty(args)
   try, g = struct(args{:});
   catch, disp('eegplugin_gibbons(): wrong syntax in function arguments'); return; end;
else
    g = [];
end;
try 
  g.called_by_temporal_distortion;
catch, 
  g.called_by_temporal_distortion=0;
end;

if g.called_by_temporal_distortion
    tempdistmenu = findobj(fig, 'tag', 'tempdist');
    uimenu( tempdistmenu, 'Label', 'Gibbons Timewarp', 'CallBack', fgibbonscmd);
else
    toolsmenu = findobj(fig, 'tag', 'tools');
    uimenu( toolsmenu, 'Label', 'Gibbons Timewarp', 'CallBack', fgibbonscmd);
end;
