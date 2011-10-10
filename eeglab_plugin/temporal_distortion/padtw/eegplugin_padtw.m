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
% $Log: eegplugin_padtw.m,v $
% Revision 1.1  2009/05/05 20:43:06  mihrke
% fixed autotools issues
%
%
%
%
function eegplugin_padtw( fig, try_strings, catch_strings, varargin);

  
% build command for menu callback
  padtwcmd = [ 'pop_padtw(ALLEEG,EEG);' ...
               'eeglab redraw;'];
  fpadtwcmd = [ try_strings.no_check padtwcmd ];
  fpadtwcmd = [fpadtwcmd  'LASTCOM = ''' padtwcmd ''';' ];
  fpadtwcmd = [fpadtwcmd  catch_strings.store_and_hist ];
  

% if called from super-plugin, put it under temporal distortion
args = varargin;
if ~isempty(args)
   try, g = struct(args{:});
   catch, disp('eegplugin_padtw(): wrong syntax in function arguments'); return; end;
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
    uimenu( tempdistmenu, 'Label', 'Hierarchical Timewarp', 'CallBack', fpadtwcmd);
else
    toolsmenu = findobj(fig, 'tag', 'tools');
    uimenu( toolsmenu, 'Label', 'Hierarchical Timewarp', 'CallBack', fpadtwcmd);
end;
