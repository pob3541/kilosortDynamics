% Create an ImecDataset pointing at a specific
% channelMapFile = 'neuropixPhase3A_kilosortChanMap.mat';
% this option works for Olaf_2018

% relevant directories
cd /media/pierreb/Biggie/NeuroPixel/2018-10-24/
cd ~/Documents/Code/neuropixel-utils/

channelMapFile = 'neuropixPhase3A_option4_kilosortChanMap.mat'; 
imec = Neuropixel.ImecDataset('/media/pierreb/Biggie/NeuroPixel/2018-10-24/Olaf_20181024_03_g0_t0.imec.ap.bin', 'channelMap', channelMapFile);


% home dir
homeDir ='/home/pierre/Code/Data/Neuropixel/2018-10-24';

channelMapFile = 'neuropixPhase3A_option4_kilosortChanMap.mat'; 
imec = Neuropixel.ImecDataset('/home/pierre/Code/Data/Neuropixel/2018-10-24/', 'channelMap', channelMapFile);




% Mark individual channels as bad based on RMS voltage
rmsBadChannels = imec.markBadChannelsByRMS('rmsRange', [3 100]);

% Inspect the raw IMEC traces
imec.inspectAP_timeWindow([200 201]); % 200-201 seconds into the recording

% Run Kilosort2
Neuropixel.runKilosort2(imec);

% Load the Kilosort2 results
ks = Neuropixel.KilosortDataset();
ks.load()

% Compute useful stats about each template and cluster
metrics = ks.getMetrics();

% Extract raw waveforms for a specific cluster id at the 24 largest amplitude channels
% Clean these waveforms by subtracting the contribution of other clusters spiking within the same time window
ss = ks.getWaveformsFromRawData('cluster_id', 255, 'num_waveforms', 100, 'best_n_channels', 24,  'subtractOtherClusters', true);


% Plot these waveforms at their physical coordinates on the neuropixel
ss.plotAtProbeLocations();


