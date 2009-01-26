function [data, times, markers] = read_raw(filename);
% read_raw(filename)
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
  
  fid = fopen(filename, 'rb');
  
  nbchan = fread(fid, 1, 'double');
  trials = fread(fid, 1, 'double');
  pnts   = fread(fid, 1, 'double');   
  num_markers= fread(fid, 1, 'double');
  times  = fread(fid, pnts, 'double');
  times = reshape( times, size(times, 2), size(times, 1));
  
  % markers
  markers = fread(fid, num_markers*trials, 'double');
  markers = reshape( markers, [num_markers, trials] );
  
  size(markers)
 
  % data
  for c=1:nbchan
      for t=1:trials
        data(c,:,t) = fread(fid, pnts, 'double');
      end;
  end;
  fclose(fid);

  return;
  
