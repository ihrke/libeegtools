function writesig(prog)

beta = 0.0;
trials = 100;
srate = 1000;
[sig real times rts urt u] = artificial_data(trials, srate, beta);
[N trials] = size(sig);


if strcmp(prog,'t_denoisestat')
    fname = sprintf('denoisestat_beta%.1f_%i_samples%i_snr+4', beta,trials, N)
    f= fopen(strcat(fname, '.dat'), 'wb');
    fwrite(f, times, 'double');
    for i = 1:trials
        fwrite(f, real(:,i), 'double');
        fwrite(f,  sig(:,i), 'double');    
    end;
    fclose(f);
elseif strcmp(prog, 't_diffmatrix')
    fname = sprintf('t_diffmatrixdata/t_diffmatrix_beta%.1f_%i_samples%i_snr-10', beta,trials, N)
    %Format of input file:
  	%	- 1st double is N
  	%	- then double n
    %   - then N*n doubles -- s_{c,i}
    f=fopen(strcat(fname, '.dat'), 'wb');
    fwrite(f, trials, 'double');
    fwrite(f, N, 'double');
    for i = 1:trials
        fwrite(f, sig(:, i), 'double');
    end;  
    fclose(f);
elseif strcmp(prog, 't_loocv')
    fname = sprintf('t_loocvdata/t_loocv_beta%.1f_%i_samples%i_snr-10', beta,trials, N)
    f=fopen(strcat(fname, '.dat'), 'wb');
    % Format of input file:
  	%	- 1st double is N
  	%	- then double n
    %   - 1st double is R_c (in real time)
    %   - following N doubles are the R_{c,i} (in real time)
    %   - then n doubles -- times array 
    %   - then n doubles -- u_c
    %   - then N*n doubles -- s_{c,i}
    zero = closest(times, 0);
    surt = closest(times,urt);
    lu = max(size(u(zero:surt)));
    fwrite(f, trials, 'double');
    fwrite(f, N, 'double');    
    fwrite(f, urt, 'double');    
    fwrite(f, rts, 'double');    
    fwrite(f, times, 'double');        
    fwrite(f, u, 'double');        
    for i = 1:trials
        fwrite(f, sig(:, i), 'double');
    end;  
    fclose(f);
elseif strcmp(prog,'t_timewarpstat')
    fname = sprintf('timewarpstat_beta%.1f_%i_samples%i_snr-10', beta,trials, N)
    f= fopen(strcat(fname, '.dat'), 'wb');
    
    zero = closest(times, 0);
    surt = closest(times,urt);
    lu = max(size(u(zero:surt)));
    fwrite(f, lu, 'double');
    fwrite(f, u(zero:surt)  , 'double');
    srts = [];

    for rt = rts
        srts = [srts closest(times, rt)];
    end;
%    srts-zero+1
    fwrite(f, srts-zero+1, 'double');
    for i = 1:trials
%        i
 %       real(zero:srts(i), i)
        fwrite(f, real(zero:srts(i), i), 'double');
    end;
    fclose(f);
elseif strcmp(prog,'t_warpavg')
    fname = sprintf('warpavg_%i_samples%i_beta%.1f', trials, N, beta)
    f= fopen(strcat(fname, '.dat'), 'wb');
    zero = closest(times, 0);
    srts = [];
    for rt = rts
       srts = [srts closest(times, rt)];
    end;

    fwrite(f, trials, 'double');
    fwrite(f, N, 'double');
    fwrite(f, zero, 'double');
    fwrite(f, srts, 'double');
    
    for i = 1:trials
        fwrite(f, real(:, i), 'double');
    end;
    fclose(f);
else
    disp('unknown prog');   
    return;
end;

save(strcat(fname, '.mat'));%, 'sig', 'real', 'times', 'srts', 'rts', 'zero', 'surt', 'urt', 'lu', 'u');
