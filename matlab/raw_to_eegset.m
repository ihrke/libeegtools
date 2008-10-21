function EEG = raw_to_eegset(filename);
% raw_to_eegset(filename);
% EEG - is a new eeglab-set
% filename - string giving name of input file
%
% NOTE:
%   * The script does not add information about electrode locations as
%     the raw-file does not know about this.
%   * The EEG-set is always epoched with markers from the raw-file.
%   * The EEG-set is not directly loaded in the GUI, you need to do
%       > EEG = raw_to_eegset( 'test.raw')
%       > [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
%       > eeglab redraw;
%     to achieve this
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
  
  EEG = eeg_emptyset;
  
  EEG.comments = sprintf('Original file: %s; created by raw_to_eegset',...
       filename);
  EEG.setname = sprintf('RAW: %s', filename);
  
  
  fid = fopen(filename, 'rb');
  
  EEG.nbchan = fread(fid, 1, 'double');
  EEG.trials = fread(fid, 1, 'double');
  EEG.pnts   = fread(fid, 1, 'double');   
  num_markers= fread(fid, 1, 'double');
  EEG.times  = fread(fid, EEG.pnts, 'double');
  reshape( EEG.times, size(EEG.times, 2), size(EEG.times, 1));
  EEG.xmin = EEG.times(1)/1000; %% eeg_checkset requiqres seconds here
  EEG.xmax = EEG.times(end)/1000;
  EEG.srate = 1000/(EEG.times(2)-EEG.times(1));
  
  % markers
  markers = fread(fid, num_markers*EEG.trials, 'double');
  
  size(markers)
 
  EEG.event = [struct()];
  for t=1:EEG.trials
    m = markers((num_markers*(t-1))+1:(num_markers*t));
    lat = m' %(EEG.times(m)')-(EEG.xmin*1000)
    
    idx = 1;
    for e=lat
        EEG.event(num_markers*(t-1)+idx).latency = (((t-1)*EEG.pnts)+e);
        EEG.event(num_markers*(t-1)+idx).type = sprintf('E%i', idx);
        EEG.event(num_markers*(t-1)+idx).code = 'Stimulus';
        EEG.event(num_markers*(t-1)+idx).duration = 0;
        EEG.event(num_markers*(t-1)+idx).channel = 0;
        EEG.event(num_markers*(t-1)+idx).time =[];
        EEG.event(num_markers*(t-1)+idx).urevent = 0;
        EEG.event(num_markers*(t-1)+idx).epoch = t;
        
        idx = idx+1;
    end;
  end;
  
  % data
  for c=1:EEG.nbchan
      for t=1:EEG.trials
          EEG.data(c,:,t) = fread(fid, EEG.pnts, 'double');
      end;
  end;
   
  fclose(fid);
  EEG = eeg_checkset( EEG, 'eventconsistency' );
  return;
  
