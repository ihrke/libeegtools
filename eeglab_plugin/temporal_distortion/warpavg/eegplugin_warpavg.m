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
% Revision 1.1  2009/01/12 11:48:00  mihrke
% warpavg
%
% Revision 1.2  2009/01/12 11:36:55  mihrke
% removed dummy
%
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
function eegplugin_warpavg( fig, try_strings, catch_strings);

%  path(path, './src')
  
% create menu
toolsmenu = findobj(fig, 'tag', 'tools');

% build command for menu callback
wavedencmd = [ 'EEG = pop_waveden(EEG); [ALLEEG EEG CURRENTSET]=eeg_store(ALLEEG,EEG);'...
        'eeglab redraw;'];
finalcmd = [ try_strings.no_check wavedencmd ];
finalcmd = [ finalcmd 'LASTCOM = ''' wavedencmd ''';' ];
finalcmd = [ finalcmd catch_strings.store_and_hist ];

warpavgcmd= [ 'EEG = pop_warpavg(EEG); [ALLEEG EEG CURRENTSET]=eeg_store(ALLEEG,EEG);'...
        'eeglab redraw;'];
fwarpavgcmd =  [ try_strings.no_check warpavgcmd ];
fwarpavgcmd = [ fwarpavgcmd 'LASTCOM = ''' warpavgcmd ''';' ];
fwarpavgcmd = [ fwarpavgcmd catch_strings.store_and_hist ];

submenu = uimenu( toolsmenu, 'label', 'Correction of Temporal Distortion');
uimenu( submenu, 'Label', 'Wavelet-Denoising', 'CallBack', finalcmd);
uimenu( submenu, 'Label', 'WarpAvg', 'CallBack', fwarpavgcmd);