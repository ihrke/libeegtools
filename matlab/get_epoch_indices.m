function idx = get_epoch_indices(EEG, eventtype);
% get the indices for the epochs where the first
% event matches one from the 'eventtype' array
  idx = [];
  for i=1:length(EEG.epoch)
    zero_idx = find([EEG.epoch(i).eventlatency{:}]==0);
    if strmatch(EEG.epoch(i).eventtype(zero_idx), eventtype, 'exact')
      idx = [idx i];      
    end;
  end;
  
