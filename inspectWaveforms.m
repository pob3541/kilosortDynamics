function [gwfparams,wF]=inspectWaveforms(clust,ksDir,sp,sr,nSpikes,qualMet)

gwfparams.dataDir = ksDir;    % KiloSort/Phy output folder
apD = dir(fullfile(ksDir, '*ap*.bin')); % AP band file from spikeGLX specifically
gwfparams.fileName = apD(1).name;         % .dat file containing the raw 
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 385;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [-40 41]; % Number of samples before and after spiketime to include in waveform

% if num wfs greater than 10k then just pull 10k
gwfparams.nWf = nSpikes;%10000; %qualMet.nSpClus(clust+1,2);                    % Number of waveforms per unit to pull out

gwfparams.spikeTimes = ceil(sp.st(sp.clu==clust)*sr); % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = sp.clu(sp.clu==clust);

% added floor to nSamp in getWaveForms

wf = getWaveForms(gwfparams);

wF.chLogic = ~isnan(qualMet.filtChs(:,clust+1));
wF.chs=qualMet.filtChs(wF.chLogic,clust+1);
wF.waveForms=squeeze(wf.waveForms);
wF.chWfMean=squeeze(mean(wF.waveForms(:,wF.chs,:),1));
wF.peak2peak=sortrows([wF.chs,max(abs(wF.chWfMean),[],2)-min(abs(wF.chWfMean),[],2)],2,'descend');
wF.bestCh=squeeze(wF.waveForms(:,wF.peak2peak(1,1),:));


end