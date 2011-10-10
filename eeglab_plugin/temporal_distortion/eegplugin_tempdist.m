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
% $Log: eegplugin_tempdist.m,v $
% Revision 1.1  2009/08/25 10:11:36  mihrke
% migrating to sf.net
%
%
%
function eegplugin_tempdist( fig, try_strings, catch_strings);
% add the paths
% -------------
  eeglabpath = which('eeglab.m');
  eeglabpath = eeglabpath(1:end-length('eeglab.m'));

  path(path, strcat(eeglabpath,'./plugins/temporal_distortion/gibbons') );
  path(path, strcat(eeglabpath,'./plugins/temporal_distortion/warpavg') );
  path(path, strcat(eeglabpath,'./plugins/temporal_distortion/padtw') );
  
  % create menu
  toolsmenu = findobj(fig, 'tag', 'tools');
  tempdistmenu = uimenu( toolsmenu, 'label', 'Correction of Temporal Distortion', ...
                         'tag', 'tempdist');
  

  % call the single plugin-fcts
  eegplugin_gibbons( fig, try_strings, catch_strings, 'called_by_temporal_distortion', 1 );
  eegplugin_padtw( fig, try_strings, catch_strings, 'called_by_temporal_distortion', 1 );
  %  eegplugin_warpavg( fig, try_strings, catch_strings, 'called_by_temporal_distortion', 1 );
