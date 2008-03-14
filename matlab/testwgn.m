noise = [1 1000];
N = 100;

s1 = rand(noise);
for i=1:N
    wc = rand(noise);
    [up sp] = ml_testwarppath(s1, wc);
    figure(2); plot(up, sp);
    pause(.1);
end;

figure(2); plot(wc);