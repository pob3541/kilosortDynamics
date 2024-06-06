function [gwfparams,wF]=inspectWaveforms(clust,ksDir,sp,sr,qualMet)
% function [gwfparams,wF]=inspectWaveforms(clust,ksDir,sp,sr,nSpikes,qualMet)

%ksDir = myKsDir;
gwfparams.dataDir = ksDir;    % KiloSort/Phy output folder
apD = dir(fullfile(ksDir, '*ap*.bin')); % AP band file from spikeGLX specifically
gwfparams.fileName = apD(1).name;         % .dat file containing the raw 
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 385;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [-40 41]; % Number of samples before and after spiketime to include in waveform

% if num wfs greater than 10k then just pull 10k
gwfparams.nWf = qualMet.nSpClus(clust+1,2);%nSpikes;%10000; %qualMet.nSpClus(clust+1,2);                    % Number of waveforms per unit to pull out

gwfparams.spikeTimes = ceil(sp.st(sp.clu==clust)*sr); % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = sp.clu(sp.clu==clust);

% added floor to nSamp in getWaveForms

%% getWaveForms
%wf = getWaveForms(gwfparams);

fileName = fullfile(gwfparams.dataDir,gwfparams.fileName);           
filenamestruct = dir(fileName);
dataTypeNBytes = numel(typecast(cast(0, gwfparams.dataType), 'uint8')); % determine number of bytes per sample
nSamp = filenamestruct.bytes/(gwfparams.nCh*dataTypeNBytes);  % Number of samples per channel
wfNSamples = length(gwfparams.wfWin(1):gwfparams.wfWin(end));
mmf = memmapfile(fileName, 'Format', {gwfparams.dataType, [gwfparams.nCh fix(nSamp)], 'x'});
channels = rmmissing(qualMet.filtChs(:,clust+1));% PB addition
chMap = channels;% readNPY(fullfile(gwfparams.dataDir, 'channel_map.npy'))+1;               % Order in which data was streamed to disk; must be 1-indexed for Matlab
nChInMap = numel(chMap);

% Read spike time-centered waveforms
unitIDs = unique(gwfparams.spikeClusters);
numUnits = size(unitIDs,1);
spikeTimeKeeps = nan(numUnits,gwfparams.nWf);
waveForms = nan(numUnits,gwfparams.nWf,nChInMap,wfNSamples);
waveFormsMean = nan(numUnits,nChInMap,wfNSamples);

for curUnitInd=1:numUnits
    curUnitID = unitIDs(curUnitInd);
    curSpikeTimes = gwfparams.spikeTimes(gwfparams.spikeClusters==curUnitID);
    curUnitnSpikes = size(curSpikeTimes,1);
    spikeTimesRP = curSpikeTimes(randperm(curUnitnSpikes));
    spikeTimeKeeps(curUnitInd,1:min([gwfparams.nWf curUnitnSpikes])) = sort(spikeTimesRP(1:min([gwfparams.nWf curUnitnSpikes])));
    for curSpikeTime = 1:min([gwfparams.nWf curUnitnSpikes])
%         tmpWf = mmf.Data.x(1:gwfparams.nCh,spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(1):spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(end));
%         waveForms(curUnitInd,curSpikeTime,:,:) = tmpWf(chMap,:);
        tmpWf = mmf.Data.x(channels,spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(1):spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(end));
        waveForms(curUnitInd,curSpikeTime,:,:) = tmpWf;
%         tmpWf = mmf.Data.x(1:gwfparams.nCh,spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(1):spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(end));
%         %tmpWf = mmf.Data.x(channels,spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(1):spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(end));
%         waveForms(curUnitInd,curSpikeTime,:,:) = tmpWf(chMap,:);
    end
    waveFormsMean(curUnitInd,:,:) = squeeze(nanmean(waveForms(curUnitInd,:,:,:),2));
    disp(['Completed ' int2str(curUnitInd) ' units of ' int2str(numUnits) '.']);
end

%%
% Package in wf struct
wf.unitIDs = unitIDs;
wf.spikeTimeKeeps = spikeTimeKeeps;
wf.waveForms = waveForms;
wf.waveFormsMean = waveFormsMean;

wF.chLogic = ~isnan(qualMet.filtChs(:,clust+1));
wF.chs=qualMet.filtChs(wF.chLogic,clust+1);
wF.waveForms=squeeze(wf.waveForms);
% wF.chWfMean=squeeze(mean(wF.waveForms(:,wF.chs,:),1));
wF.chWfMean=squeeze(mean(wF.waveForms,1));
wF.peak2peak=sortrows([wF.chs,max(abs(wF.chWfMean),[],2)-min(abs(wF.chWfMean),[],2)],2,'descend');
% wF.bestCh=squeeze(wF.waveForms(:,wF.peak2peak(1,1),:));
bestCh=wF.peak2peak(1,1);
wF.bestCh=squeeze(wF.waveForms(:,find(bestCh==chMap),:));



% wF.chWfMean=squeeze(mean(wF.waveForms(:,wF.chs,:),1));
% % consider getting rid of the absolute value for the min as likely changes
% % it to zero for some - like in the matlab function for it
% disp('consider changing to not the absolute value in this case')
% wF.peak2peak=sortrows([wF.chs,max(abs(wF.chWfMean),[],2)-min(abs(wF.chWfMean),[],2)],2,'descend');
% wF.bestCh=squeeze(wF.waveForms(:,wF.peak2peak(1,1),:));

wF.chWfMean=squeeze(wf.waveFormsMean);
% wF.chWfMean=squeeze(mean(wF.waveForms(:,wF.chs,:),1));
% wF.peak2peak=sortrows([wF.chs,max(abs(wF.chWfMean),[],2)-min(abs(wF.chWfMean),[],2)],2,'descend');
wF.peak2peak=sortrows([wF.chs,max(wF.chWfMean,[],2)-min(wF.chWfMean,[],2)],2,'descend');
wF.bestCh=squeeze(wF.waveForms(:,9,:));%wF.peak2peak(1,1),:));



end