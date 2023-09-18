% run kilosort

cd /home/pierre/Code/Kilosort
kilosort

% set advanced options in GUI creates ks object

%cd /home/pierre/Code/Kilosort
kilosort

% run kilosort automatically 
dataFile = '';
dPath = '/home/pierre/Code/Data/Neuropixel/';
oPath = '';
probeLayout = '';
nChs = 384; % Kilosort magically guesses a number
sf = 30000;
tRng = [0 Inf];
nBlocksRegistration = 5; % ??
Threshold = [9 9]; % ??
Lambda = 10; %??
AUCsplits=[0.9]; %??

ks.ops.chanMap
ks.ops.fbinary=[dPath,'2018-10-24/Olaf_20181024_02_g0_t0.imec.ap.bin'];
