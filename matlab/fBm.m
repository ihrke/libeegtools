function [ noise ] = fBm( beta, n )
%FBM generate a n-vector of 1/f^beta noise


% 	real[0] = 0;
% 	imag[0] = 0;
% 	for (i=1;i<=N/2;i++) {
% 		mag = pow(i+1.0,-beta/2) * RandomGaussian(0.0,1.0); // Note to self a number of years later, why "i+1"
% 		pha = TWOPI * RandomUniform();
% 		real[i] = mag * cos(pha);
% 		imag[i] = mag * sin(pha);
% 		real[N-i] =  real[i];
% 		imag[N-i] = -imag[i];
% 	}
% 	imag[N/2] = 0;
% 
% 	FFT(-1,TWOPOWER,real,imag);

re = zeros(n,1);
im = zeros(n,1);
for i=1:n/2
    mag = (i+1)^(-beta/2) * randn(1);
    pha = 2*pi * rand(1);
    re(i) = mag * cos(pha);
    im(i) = mag * sin(pha);
    re(n-i) = re(i);
    im(n-i) = im(i);
end;
im(round(n/2)) = 0;


noise = real(ifft(complex(re, im)));
return;