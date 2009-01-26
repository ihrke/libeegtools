function [real_erp single_trials_nonoise] = write_art_data_raw(filename, trials);
plotit =0;  

[single_trials_noisy single_trials_nonoise range rts real_rt real_erp] =  artificial_data(trials, ...
                                                  500, 0);
rts

d = single_trials_nonoise;

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
fwrite(fid, range, 'double');
  
% markers in sampling points
for t=1:trials
    % stimulus onset
    fwrite(fid, closest(range, 0), 'double');
    % reaction
    fwrite(fid, closest(range, rts(t)), 'double');
    closest(range, rts(t))
end;

% data
for i=1:trials
  fwrite(fid, d(:,i), 'double');
end;
fclose(fid);

if plotit
  figure;
  subplot(2,2,1);
  plot(range, erp1, 'k', 'LineWidth', 4);
%  subplot(2,2,3);
%  plot(range, erp2, 'k', 'LineWidth', 4);

  subplot(2,2,[2 4]);
  plot(range, csnew);  
end;