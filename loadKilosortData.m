% loadKilosortdata

myKsDir = '/home/pierre/Code/Data/NeuroPixel/2018-10-24/output';
cd(myKsDir)
sp = loadKSdir(myKsDir); % from spikes repository

% crashes matlab

% ds='Olaf_20181024_02_g0_t0.imec.ap.bin';
% fileID = fopen(ds);
% A = fread(fileID);

% load in waveforms


% similarity scores to look into ....
load rez.mat
sort(rez.simScore(:),'descend')
simScoreLogic=rez.simScore>0.1 & rez.simScore<0.9999999;
idx=find(simScoreLogic);
scores=rez.simScore(simScoreLogic);




sp.st % spike times in seconds

sp.cgs % manual labels assigned to the units

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


%%

spikeTimes=sp.st;
spikeTemps = sp.spikeTemplates;
clusters=sp.clu;
ycoords = sp.ycoords;


% sum number of instances of each cluster
for i =0 : max(clusters)
storClu(:,i+1)=clusters==i;
end
nNans=max(sum(storClu));

% spike times organized by cluster
for i = 1: size(storClu,2)
tmpSpClus=spikeTimes(storClu(:,i));
nanFill=nNans-length(tmpSpClus);
spClus(:,i)=[tmpSpClus;nan(nanFill,1)];
end

% ISI calculation
sr= sp.sample_rate;
for i = 1:10%size(spClus,2)
tmpSpClus=spClus(:,i);
ISI=diff(tmpSpClus(~isnan(spClus(:,i))))*1000; % in ms
figure;
%boxplot(ISI);
histogram(ISI(ISI<1000),'Normalization','probability')
boxchart(ISI,'MarkerStyle','none')



%histogram(ISI); xlabel('Time (ms)');
%ylabel ('ISI count')

% nan logic to store all the ISIs - could do a boxplot
end


% the rule of thumb is that if of the ISIs less than 1.5% violate the ISI
% (less  than 1.5 ms) then it is labelled as a single unit



% how to load in the raw data 
% plot a raw data channel trace - just load the binary data


% which channels are these clusters on 
j=1;
for i =1:size(exTemp,2)
nTPs=size(exTemp,1);

   % can use any here to find the zero vectors
    if ~isempty(find(exTemp(:,i), 1))
        store(:,j)=exTemp(:,i);
        j=j+1;
    end

end

for i =1:size(store,2)
figure;plot(store(:,i))
end

%% play around with the structure
% could potentially calculate the pc features myself

% sp.temps
% how are they calculated
for i =40:50
figure;plot(sp.temps(10,:,i))
end
exTemp=squeeze(sp.temps(1,:,:));



% sp.spikeTemplates - at the spike time which template it is
sp.spikeTemplates
