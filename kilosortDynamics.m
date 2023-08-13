%What data quality metrics Can we use?

% nSpikes count spikes from a single session nSpikes/sess
% 

% What is the model of our probe? 

% install MatNWB  
addpath(genpath(pwd));
generateCore();


% load dataset in 
path='/home/pierre/Code/Data/kilosortData/';
cd(path)
fileName='Steinmetz2019_Theiler_2017-10-11.nwb';
 

% read NWB file
read_nwbfile = nwbRead([path,fileName]);

% load ChanMap.mat code
addpath map_files/
labProbe = 'neuropixPhase3A';
load([labProbe, '_kilosortChanMap.mat']);


imec = Neuropixel.ImecDataset('/data/raw_datasets/neuropixel_01_g0_t0.imec.ap.bin');

