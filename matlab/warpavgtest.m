[s us r rts rt u] = artificial_data(100, 1000, 0);

w=ml_warpavg(r, s, rts);
figure; 
plot(r, mean(s,2), 'y');
hold on;
plot(r, mean(us,2), 'b');
plot(r, w,'k');
plot(r, u, 'r');