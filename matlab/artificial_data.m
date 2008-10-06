function [single_trials_noisy single_trials range rts real_rt erp trans] = ...
    artificial_data(trials, srate, beta);
% function [single_trials_noisy single_trials range rts real_rt erp] = ...
%    artificial_data(trials, srate, beta);
%clear all;
%close all;

real_rt = 800; % "real" (u) reaction time
sd_rt = 100;
%srate = 1000; % sampling rate (Hz)

% sets: 0.8/1000; 0/100; 2/5000 
noisebeta= 0;  
%noiseamp = 100;
%trials = 100;
plotit = 0;

rts = [];
rand('state',sum(100*clock))
step = 1000/srate;
if plotit,figure; end;
range = -500:step:1998; % in ms

% ERP - Spline interpolation
%fixp = [0 0; 100 10; 200 -10; 300 20; 500 0; 800 0;]
fixp1 = [
    -500 0.8;
    -450 -0.5;
    -380 0.4;
    -333 0.2;
    -150 2;
    0 0;
    50 -3;
    100 4.3;
    200 -13;
    250 -2;
    300 -3;
    500 15;
    600 8;
    700 10;
    900 -2;
    1000 0.5;
    1200 -3;
    1500 2;
    1950 -2;
    2000 0;
];

fixp2 = [
    -500 0.8;
    -450 -0.5;
    -380 0.4;
    -333 0.2;
    -150 2;
    0 0;
    50 3;
    100 -4.3;
    200 13;
    250 2;
    300 3;
    500 -15;
    600 8;
    700 10;
    900 -2;
    1000 0.5;
    1200 -3;
    1500 2;
    1950 -2;
    2000 0;
];

fixp=fixp1;
erp = spline(fixp(:,1), fixp(:,2), range);
noiseamp = max(erp)*100;

if plotit
    subplot(3,2,1);
    plot(range, erp);
    title('Simulated "real" ERP');
    xlim([range(1) range(end)]);
    subplot(3,2,3);
end;

% generating the phi_i

% points in sampling notation
zero = closest(range, 0);
rresp = closest(range, real_rt);


single_trials = [];
for i = 1:trials
    rt = real_rt + sd_rt*randn(1); % gaussian rt
    rts = [rts rt];
    k = 1;
    trans = interp1([range(1) 0 real_rt range(end)],...
        [range(1) 0 rt range(end)],...
        range);
     trans(zero:rresp-1) = get_monotonous_function([0 0],...
          [real_rt rt], step);
        
    terp = interp1(trans, erp, range)';
    alpha = rand+0.5;
%    alpha=abs(randn+1);
    terp = alpha*terp;
    if plotit
        plot(range, terp, 'r');
        if i==1
            hold on
        end;
    end;
    single_trials = [single_trials terp];
end;
[N p] = size(single_trials);
single_trials_noisy=[];
for i=1:p
    % 1150 for beta=0
    % 5000 for beta=0.5
    % 15000 for beta=1
    % 32000 for beta=1.5
    % 55000 for beta=2.0
    if beta==0
        namp=1000;
    elseif beta==0.5
        namp=5000;
    elseif beta==1.0
        namp=15000;
    elseif beta==1.5
        namp=32000;
    elseif beta==2.0
        namp=55000;
    else
        disp('ERROR');
        return;
    end;
    noise = (fBm(beta, N)).*namp;
    single_trials_noisy(:,i) = single_trials(:,i)+noise;%30*randn(N,1);%noise;
end;

if plotit
    plot(range, trans, 'b');
    hold on;
    plot(trans, range, 'r');
    xlim([range(1) range(end)]);
    hold off;
end;

if plotit
    title(sprintf('Single trials, generated with gaussian rts with sd=%i', sd_rt));
    xlim([range(1) range(end)]);
    hold off;

    subplot(3,2,2);
    plot(range, trans, 'b');
    hold on;
    plot(trans, range, 'r');
    title(sprintf('Sample phi_i and its inverse for rt=%3.2f', rt));
    xlim([range(1) range(end)]);
    hold off;

    subplot(3,2,4);
    plot(range, mean(single_trials, 2));
    hold on
    plot(range, mean(single_trials,2)-erp', 'k');
    title('Mean of single trials (black curve is error)');
    xlim([range(1) range(end)]);
    hold off;


    subplot(3, 2, 5);
    [N p] = size(single_trials);
    plot(range, single_trials_noisy(:,i));
    title(sprintf('Sample Trial with 1/f noise added, beta=%.1f',noisebeta));
    xlim([range(1) range(end)]);
    hold off;

    subplot(3, 2, 6);
    [N p] = size(single_trials);
    plot(range, mean(single_trials_noisy, 2));
    hold on;
    %plot(range, mean(single_trials_noisy,2)-erp', 'k');
    title(sprintf('Mean of trials with 1/f noise added, beta=%.1f',noisebeta));
    xlim([range(1) range(end)]);
    hold off;
end;
