N=20;
srate=500;

[s u r rt R] = artificial_data(N, srate, 0);
zero = closest(r, 0)
rt
srt=[];
for t=rt
  srt = [srt closest(r, t)];
end;
srt

u;
w = ml_padtw(u, zero, srt);

figure; 
for i=1:N
  plot(r, u(:,i));
  if i==1, hold on; end;
end;
plot(r, w, 'r');