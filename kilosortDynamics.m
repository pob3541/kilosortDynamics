%What data quality metrics Can we use?

% nSpikes count spikes from a single session nSpikes/sess
% 

% What is the model of our probe? 

% install MatNWB  
addpath(genpath(pwd));
generateCore();


% load datasets in 
% steinmetz-et-al-2019
path='/home/pierre/Code/Data/kilosortData/';
cd(path)
fileName='Steinmetz2019_Theiler_2017-10-11.nwb';

% read NWB file
nwbfile = nwbRead([path,fileName]);
units=nwbfile.units;
spikeTimes=units.spike_times;
spikeTimes.data

% paulk-et-al-2023
dataPath ='/home/pierre/Code/Data/paulk-et-al-2023/pt2/Pt02.imec0.ap.bin';
chPath ='/home/pierre/Code/Data/paulk-et-al-2023/chMap';
cd(path)
%paulkProbe ='chanMap_DREDgeAligned_Pt02.mat';

% load ChanMap.mat code
addpath /home/pierre/Code/neuropixel-utils/map_files/
guessPaulkProbe = 'neuropixPhase3B1';
load([guessPaulkProbe, '_kilosortChanMap.mat']);


%construct neuropixel imec dataset by pointing to path
channelMapFile = guessPaulkProbe ;
channelMapFile = 'neuropixPhase3A_kilosortChanMap.mat';

imec = Neuropixel.ImecDataset('/home/pierre/Code/Data/paulk-et-al-2023/pt2/Pt02.imec0.ap.bin', 'channelMap', channelMapFile);

% Create an ImecDataset pointing at a specific
imec = Neuropixel.ImecDataset('/data/raw_datasets/neuropixel_01.imec.ap.bin', 'channelMap', channelMapFile);
load





















%imec = Neuropixel.ImecDataset('/data/raw_datasets/neuropixel_01_g0_t0.imec.ap.bin');

