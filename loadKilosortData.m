%% loadKilosortdata

myHomeKsDir = '/home/pierre/Code/Data/NeuroPixel/2018-10-24/output';
myWorkKsDir = '/media/pierreb/Biggie/NeuroPixel/2018-10-24/';

%workCodeDir 
cd '~/Documents/Code/kilosortDynamics'

cd(myWorkKsDir)
sp = loadKSdir(myWorkKsDir); % from spikes repository

%% calculate quality metrics - to determine which units waveforms to plot

[spTimes,clusIdx,spTemp,qualMet]=calcISI(sp);

table(qualMet.nSpClus)
% which channels are these clusters on 
nClusters=qualMet.nClusts;
for clus =1: nClusters
%clus=1;
exTemp=squeeze(sp.temps(clus,:,:));

j=1;
for i =1:size(exTemp,2)
nTPs=size(exTemp,1);
chLogic(i,1)=0;

   % can use any here to find the zero vectors
    if ~isempty(find(exTemp(:,i), 1))
        nonZeroTemp(:,j)=exTemp(:,i);
        chLogic(i,1)=1;
        j=j+1;
    end

end
channels =1:383;
filtChs=channels(logical(chLogic))';



% will plot the template from certain channels
figure;
hold on;
for i =1:length(filtChs)
plot(sp.temps(clus,:,filtChs(i)))
end
end

%%
% load in waveforms 
% added floor to nSamp in getWaveForms
myKsDir=myWorkKsDir;
clust=66;

gwfparams.dataDir = myKsDir;    % KiloSort/Phy output folder
apD = dir(fullfile(myKsDir, '*ap*.bin')); % AP band file from spikeGLX specifically
gwfparams.fileName = apD(1).name;         % .dat file containing the raw 
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 385;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
gwfparams.spikeTimes = ceil(sp.st(sp.clu==clust)*30000); % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = sp.clu(sp.clu==clust);

wf = getWaveForms(gwfparams);

% the cluster of spikes seems to move significantly over time
figure; 
imagesc(squeeze(wf.waveFormsMean))
set(gca, 'YDir', 'normal'); xlabel('time (samples)'); ylabel('channel number'); 
colormap(colormap_BlueWhiteRed); caxis([-1 1]*max(abs(caxis()))/2); box off;

%plot individual waveforms
waveForms=squeeze(wf.waveForms);

figure;
chs =28;
nSpikes=81;
for i= 1:nSpikes 
hold on;
plot(waveForms(:,chs,i))
end

%%







% how to load in the raw data 
% plot a raw data channel trace - just load the binary data





% play around with the structure
% could potentially calculate the pc features myself




% sp.spikeTemplates - at the spike time which template it is
sp.spikeTemplates

% similarity scores to look into ....
load rez.mat
sort(rez.simScore(:),'descend')
simScoreLogic=rez.simScore>0.1 & rez.simScore<0.9999999;
idx=find(simScoreLogic);
scores=rez.simScore(simScoreLogic);

%ks.rez are the results from kilosort

%results xc, yc are coordinates of each channel on the probe




sp.cgs % manual labels assigned to the units






%% need the pc features to use these functions here
% also Dan O'Shea also has tools to use here

% trying to the driftmap and the output but missing pc_features making the
% function not useable
loadPCs =false;
empt={1};
[spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(myKsDir);



temps = sp.temps;
winv = sp.winv;
tempScalingAmps = sp.tempScalingAmps;

[spikeAmps, ~, templateYpos, tempAmps, tempsUnW] = ...
  templatePositionsAmplitudes(temps, winv, ycoords, spikeTemps, tempScalingAmps);

depthBins = 0:40:3840;
ampBins = 0:30:min(max(spikeAmps),800);
recordingDur = sp.st(end);

% need pc features again for spikeDepths

[pdfs, cdfs] = computeWFampsOverDepth(spikeAmps, spikeDepths, ampBins, depthBins, recordingDur);
plotWFampCDFs(pdfs, cdfs, ampBins, depthBins);
