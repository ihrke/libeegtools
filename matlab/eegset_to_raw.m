function eegset_to_raw(EEG, num_markers, events, filename);
%  eegset_to_raw(EEG, num_markers, events, filename);
% EEG - is an eeglab-set
% num_markers - how many markers per trial?
% events - is an array giving the indices from the
%          EEG.epoch field to be used
% filename - string giving name of output file
%
% Output file-format:
% (each segment in the table is a 64 bits double):
%	  1) number of channels
%	  2) number of trials
%	  3) number of samples per trial segment
%	  4) number of markers per trial
%	  5) n-doubles giving the times-array (for each sample from 1...n this array
%	     gives the corresponding time in ms)
% 	  6) num_markers*num_trials-doubles giving markers in sampling-point
% 	  units [1,...,n]
%	  7) raw EEG-data in the format of:
%	      - channels x trial x samples
%			- i.e. first all trials of the first channel one after the other, than
%			  the 2nd and so on
  fid = fopen(filename, 'wb');
  fwrite(fid, EEG.nbchan, 'double');
  fwrite(fid, max(size(events)), 'double');
  fwrite(fid, EEG.pnts, 'double');
  fwrite(fid, num_markers, 'double');
  fwrite(fid, EEG.times, 'double');
  
  % markers
  ms = EEG.epoch(events);
  for i=1:max(size(events))
    for j=1:num_markers
      t = ms(i).eventlatency(j);  
      t = closest(EEG.times, t{:})
      fwrite(fid, t, 'double');
    end;
  end;
  
  % data
  for c=1:EEG.nbchan
    for i=events
      fwrite(fid, EEG.data(c,:,i), 'double');
    end;
  end;
  
  fclose(fid);
  
