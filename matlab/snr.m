function s = snr( realsig, sig )
%SNR computes signal to noise-ratio (Wang et al., formula (9))

s = 10*log10(sum(realsig.^2) / sum((realsig-sig).^2));

return;