%% you need to change most of the paths in this block

addpath(genpath('/home/pierreb/Documents/Code/Kilosort')) % path to kilosort folder
addpath('/home/pierreb/Documents/Code/npy-matlab') % for converting to Phy
rootZ = '/media/pierreb/Biggie/2023-09-18_09-31-11/Record Node 101/experiment1/recording1/continuous/Acquisition_Board-100.Rhythm Data'; % the raw data binary file is in this folder'
rootH = '/media/pierreb/Biggie/tempBin'; % path to temporary binary file (same size as data, should be on fast SSD)
output='/home/pierreb/Documents/Code/sortOutput';
pathToYourConfigFile = '/home/pierreb/Documents/Code/kilosortDynamics'; % take from Github folder and put it somewhere else (together with the master_file)
chanMapFile = 'plexon32ChanMap.mat';

ops.trange    = [0 Inf]; % time range to sort
ops.NchanTOT  = 40; % 32-40 including the electrodes and ADC channels% total number of channels in your recording

run(fullfile(pathToYourConfigFile, 'plexonConfig.m'))
ops.fproc   = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD
ops.chanMap = fullfile(pathToYourConfigFile, chanMapFile);
%% this block runs all the steps of the algorithm
% main parameter changes from Kilosort2 to v2.5
ops.sig        = 20;  % spatial smoothness constant for registration
ops.fshigh     = 300; % high-pass more aggresively
ops.nblocks    = 5; % blocks for registration. 0 turns it off, 1 does rigid registration. Replaces "datashift" option. 

% main parameter changes from Kilosort2.5 to v3.0
ops.Th       = [9 9];

% is there a channel map file in this folder?
% fs = dir(fullfile(rootZ, 'chan*.mat'));
% if ~isempty(fs)
%     ops.chanMap = fullfile(rootZ, fs(1).name);
% end

% find the binary file
fs          = [dir(fullfile(rootZ, '*.bin')) dir(fullfile(rootZ, '*.dat'))];
ops.fbinary = fullfile(rootZ, fs(1).name);

rez                = preprocessDataSub(ops);

%% try to solve this or figure out if true is necessary
rez                = datashift2(rez, 0); % changed do_correction from 1 to 0

[rez, st3, tF]     = extract_spikes(rez);

rez                = template_learning(rez, tF, st3);

[rez, st3, tF]     = trackAndSort(rez);

rez                = final_clustering(rez, tF, st3);

rez                = find_merges(rez, 1);

rootZ = fullfile(rootZ, 'kilosort3');
mkdir(rootZ)
rezToPhy2(rez, rootZ);

% store more accessibly
%rezToPhy2(rez,output)

%% 
