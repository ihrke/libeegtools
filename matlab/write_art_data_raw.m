function [realclust cs csnew map] = write_art_data_raw(filename);
trials = 100;  
plotit =1;  

[single_trials_noisy c1 range rts real_rt erp1] =  artificial_data(trials, ...
                                                  500, 0);
rts
%[single_trials_noisy c2 range rts real_rt erp2] =  artificial_data2(trials/2, ...
%                                                  200, 0);

map = randperm(trials);
cs = [c1];

%alpha = 2*rand(trials, 1);
for i=1:trials
  cs(:,i)=cs(:,i);%*alpha(i);
end;

csnew = cs(:,map);
%realclust = (~(map<=trials/2))+1;

% num channels
fid = fopen(filename, 'wb');
fwrite(fid, 1, 'double');

% num trials
fwrite(fid, trials, 'double');

% num samples per trial
fwrite(fid, size(cs, 1) , 'double');

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
  fwrite(fid, csnew(:,i), 'double');
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