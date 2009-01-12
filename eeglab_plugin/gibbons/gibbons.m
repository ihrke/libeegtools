% gibbons() -  Calculates Response-Time-Corrected Average as 
%                   detailed in Gibbons+Stahl (2007) 'Response-time
%                   Corrected averaging of event-related potentials'
%                   in Clinical Neurophysiology, 118, 197-208
%
% Usage:
%   >>  geeg = gibbons( eeg, times, trialRT, k );
%
% Inputs:
%   eeg     - to-be-corrected data in the format:
%             Electrode x Time x trial
%             such that eeg(3, 400, 4) is the 3rd electrode at sampling point
%             400 of trial 4
%   times   - array of original times
%   trialRT - reaction time fo each trial (vector)
%   k       - parameter for the warping (see Gibbons+Stahl 2007)
%    
% Outputs:
%   ceeg    - corrected EEG data (same format as input)
%
% See also: 
%   POP_GIBBONS, EEGLAB 

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

function ceeg = gibbons( eeg, times, trialRT, k )

if nargin < 2
	help gibbons;
	return;
end;

meanRT = mean(trialRT);
zero_idx = closest(times, 0);
ntimes = times(1:zero_idx);

[chans samples trials]=size(eeg);


for i=1:trials
  trialRT_idx=closest( times, trialRT(i) );
  t = ntimes;
  rt = times( zero_idx:trialRT_idx );
 
  t(zero_idx:trialRT_idx) = (rt) + ((rt).^k)./(trialRT(i).^k) .* ...
      (meanRT-trialRT(i));
  
  t(trialRT_idx:samples)=1;
  t(trialRT_idx:samples)=linspace(meanRT, times(end),...
      size(t(trialRT_idx:samples),2));
  for c=1:chans
      ceeg(c,1:samples,i) = interp1(t, eeg(c,:,i), times, 'linear');
  end;
end;



return;