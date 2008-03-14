function r = rmse( realsig, sig )
%SNR computes signal to noise-ratio (Wang et al., formula (9))

r = sqrt( 1/length(realsig) * sum( (realsig-sig).^2));

return;
