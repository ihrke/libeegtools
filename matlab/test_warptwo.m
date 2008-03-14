[si ui times Ri R u] = artificial_data(100, 500, 0);
scale = 2.0;

sRi1 = closest(times, Ri(1));
sRi2 = closest(times, Ri(2));
zero = closest(times, 0);
s1 = ui(:,1)';
s2 = 2*ui(:,2)';
grey = [0.9 0.9 0.9];
%s2 = scale*s2;

wc = ml_warptwo(s1, s2, [zero sRi1 sRi2]);
%[up sp] = ml_testwarppath(s1, s2);

s1p = s1(zero:sRi1);
s2p = s2(zero:sRi2);
[up sp] = ml_testwarppath(s1p, s2p);

figure(1);
plot(times, s1, 'Color', grey, 'Linewidth', 2);
hold on;
plot([Ri(1) Ri(1)], [-40 40], '--');
plot(times, s2, 'Color', grey, 'Linewidth', 2);
plot([Ri(2) Ri(2)], [-40 40], '--');
plot(times, wc, 'r', 'Linewidth', 4);
hold off;

figure(2)
plot(up(1:end-2), sp(1:end-2));


%plot(3);
%plot(up2, sp2);

%plot(4); 
%plot(times(zero:sRi1), s1p, 'b');
%hold on
%plot(times(zero:sRi2), s2p, 'g');
%plot();

hold off;