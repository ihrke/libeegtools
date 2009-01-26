function plot_spectrum( v, sampfreq, color )

L = max(size(v));
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(v,NFFT)/L;
f = sampfreq/2*linspace(0,1,NFFT/2);

% Plot single-sided amplitude spectrum.
plot(f,log(abs(Y(1:NFFT/2)).^2), color) 
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
