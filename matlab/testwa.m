trials = 100;
srate = 1000;
beta = 1.0;

[s us r rts rt u] = artificial_data(trials, srate, beta);
wa = ml_warpavg(r, s, rts, 6, 1, 1, 2, 1, 1, 100, 100, 10, 5);
%                          L, lambda, eta, sige, theta1/2, tau, b, m, cc
for i=1:trials
    snew(i,:) = ml_denoise(s(:,i), 6, 1, 1, 2);
end;
snr(u, wa)
snr(u, mean(s,2)')
snr(u, mean(us,2)')
snr(u, mean(snew,1))
figure;
plot(r, mean(s,2), 'y'); hold on;
plot(r, u, 'r'); 
plot(r, wa, 'k');
plot(r, mean(snew,1), 'c');
plot(r, mean(us,2),'b');