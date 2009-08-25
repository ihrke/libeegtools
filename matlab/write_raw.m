function  write_raw(data, times, markers, filename);
% write_raw(filename)
% filename - string giving name of input file
%
% NOTE:
%   * 
%
% Input file-format (binary):
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
  
  [nchan n ntrials] = size(data)
  [num_markers tmp] = size(markers)
  
  fwrite(fid, nchan, 'double');
  fwrite(fid, ntrials, 'double');
  fwrite(fid, n, 'double');  
  fwrite(fid, num_markers, 'double');

  % times
  fwrite(fid, times, 'double');

  % markers
  fwrite(fid, markers(:), 'double');
  
  % data
  for c=1:nchan
      for t=1:ntrials
          fwrite(fid, data(c,:,t), 'double');
      end;
  end;
  fclose(fid);

  return;
  
