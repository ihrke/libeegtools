function [realclust cs csnew map] = write_art_clusters(filename);
trials = 20;  
plotit =1;  

[single_trials_noisy c1 range rts real_rt erp1] =  artificial_data(trials/2, ...
                                                  200, 0);

[single_trials_noisy c2 range rts real_rt erp2] =  artificial_data2(trials/2, ...
                                                  200, 0);

map = randperm(trials);
cs = [c1 c2];

alpha = 2*rand(trials, 1);
for i=1:trials
  cs(:,i)=cs(:,i)*alpha(i);
end;

csnew = cs(:,map);
realclust = (~(map<=trials/2))+1;

fid = fopen(filename, 'wb');
fwrite(fid, 1, 'double');
fwrite(fid, trials, 'double');
fwrite(fid, size(cs, 1) , 'double');
fwrite(fid, 0, 'double');
fwrite(fid, range, 'double');
  
  
% data
for i=1:trials
  fwrite(fid, csnew(:,i), 'double');
end;
fclose(fid);

if plotit
  figure;
  subplot(2,2,1);
  plot(range, erp1, 'k', 'LineWidth', 4);
  subplot(2,2,3);
  plot(range, erp2, 'k', 'LineWidth', 4);

  subplot(2,2,[2 4]);
  plot(range, csnew);  
end;