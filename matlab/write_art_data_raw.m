function [real_erp single_trials_nonoise single_trials_noisy times rts] = write_art_data_raw(filename, trials, srate);
plotit =0;  

[single_trials_noisy single_trials_nonoise times rts real_rt real_erp] =  artificial_data(trials, ...
        srate, 0);
mean(rts)
real_rt


d = single_trials_nonoise;
%d = single_trials_noisy;

alpha = (2*rand(trials, 1));
for i=1:trials
  d(:,i)=d(:,i)*alpha(i);
end;

% num channels
fid = fopen(filename, 'wb');
fwrite(fid, 1, 'double');

% num trials
fwrite(fid, trials, 'double');

% num samples per trial
fwrite(fid, size(d, 1) , 'double');

% nummarkers/trial
fwrite(fid, 2, 'double');

%times-array
fwrite(fid, times, 'double');
  
% markers in sampling points
for t=1:trials
    % stimulus onset
    fwrite(fid, closest(times, 0), 'double');
    % reaction
    fwrite(fid, closest(times, rts(t)), 'double');
    closest(times, rts(t));
end;

% data
for i=1:trials
  fwrite(fid, d(:,i), 'double');
end;
fclose(fid);

if plotit
  figure;
  subplot(2,2,1);
  plot(times, erp1, 'k', 'LineWidth', 4);
%  subplot(2,2,3);
%  plot(times, erp2, 'k', 'LineWidth', 4);

  subplot(2,2,[2 4]);
  plot(times, csnew);  
end;