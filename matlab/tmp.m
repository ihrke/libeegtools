X(1,:,:) = [ 1 2; 3 4; 5 6; 7 8; 9 10 ];   
X(2,:,:) = [ 11 12; 13 14; 15 16; 17 18; 19 20];

markers = [ 0 2 3; 
            1 3 1];

% ml_padtw( X, markers );
% 
% eeglab;
% EEG = pop_loadset( 'filename', 'eeg061206_1_resampled500hz_filtered_CO.set', 'filepath', '/scratch.local/ihrke/data/voicekey/eeglab/final/');
% EEG = eeg_checkset( EEG );
% eeglab redraw;

% settings = struct();
% settings.sampling_rate=500;
% settings.corner_freqs=[0 20];
% %settings.regularize='none';
% settings.pointdistance='euclidean_derivative';
% [avg, m] = ml_padtw( double(EEG.data(1,:,:)), [], settings )
% 


function d=ldist( Q1, Q2, P )
    d = abs(det([Q2-Q1;P-Q1]))/norm(Q2-Q1);
