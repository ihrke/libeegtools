function eegset_to_raw(EEG, num_markers, trials, filename);
%  eegset_to_raw(EEG, num_markers, events, filename);
% EEG - is an eeglab-set
% num_markers - how many markers per trial?
% trials - is an array giving the indices from the
%          EEG.epoch field to be used; [] - all trials
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
 
  if isempty(trials)
    trials = 1:length(EEG.epoch);
  end;
  
  fid = fopen(filename, 'wb');
  fwrite(fid, EEG.nbchan, 'double');
  fwrite(fid, max(size(trials)), 'double');
  fwrite(fid, EEG.pnts, 'double');
  fwrite(fid, num_markers, 'double');
  fwrite(fid, EEG.times, 'double');
 
  
  % markers
  ms = EEG.epoch(trials);
  for i=1:max(size(trials))
    for j=1:num_markers
      j
      t = ms(i).eventlatency(j)
      t = closest(EEG.times, t{:})
      fwrite(fid, t-1, 'double');
    end;
  end;
  
  % data
  for c=1:EEG.nbchan
    for i=trials
      fwrite(fid, EEG.data(c,:,i), 'double');
    end;
  end;
  
  fclose(fid);
  
  fid=fopen([filename(1:length(filename)-3) 'info'], 'w');
  fwrite(fid, sprintf('Generated from EEGlab-set:\n  %s\n', EEG.setname));
  fwrite(fid, sprintf('File: %s\n', EEG.filename));   
  fwrite(fid, sprintf('Path: %s\n', EEG.filepath));
  fwrite(fid, sprintf('Comments: %s\n', char(EEG.comments)));
  
