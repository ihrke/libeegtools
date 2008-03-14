clear all; close all;
N= 100;
[si ui times Ri R u] = artificial_data(N, 500, 0);

for i=1:N
    sRi(i) = closest(times, Ri(i));
end;
zero = closest(times, 0);

markers = [repmat(zero,[N 1]) sRi'];

uis = ml_denoise(si, 6, 0);
wa = ml_warpavg(uis, markers);


figure(2)
for i=1:N
    plot(times, si(:,i), 'Color', [0.9 0.9 0.9]);
    if i==1
        hold on;
    end;
end;
plot(times, u, 'r', 'LineWidth', 3);
plot(times, wa, 'b', 'LineWidth', 2);

plot(times, mean(uis, 2)', 'k');
hold off;
    
snrwa = snr(u, wa)
snrmean=snr(u, mean(uis, 2)')
rmsewa = rmse(u, wa)
rmsemean=rmse(u, mean(uis, 2)')

% 
% figure(1);
% plot(times, s1, 'b');
% hold on;
% plot([Ri(1) Ri(1)], [-10 10], 'b');
% plot(times, s2, 'm');
% plot([Ri(2) Ri(2)], [-10 10], 'm');
% plot(times, wc, 'r');
% hold off;
