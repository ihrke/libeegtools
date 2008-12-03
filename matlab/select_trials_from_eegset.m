function NEEG = select_trials_from_eegset(EEG, trials);
%  eegset
% EEG - is an eeglab-set
% trials is either
%      1) a numeric array with the numbers of trials to use
%      2) a string-array giving the event-code used for clustering
%
  ctrials=[];
  if iscellstr(trials)
    disp('Trials is a cell-string variable, looking for events...');
    for i=1:length(EEG.epoch)
      stimI = find([EEG.epoch(i).eventlatency{:}]==0);
      if strmatch(EEG.epoch(i).eventtype(stimI), trials, 'exact')
        ctrials = [ctrials i];
      end;                                                                                                                                     
    end;
  elseif isnumeric(trials)
      ctrials = trials;
  else 
    disp('ERROR, wrong trial variable');
    return;
  end;
  
  NEEG = pop_select(EEG, 'trial', ctrials);
  return;