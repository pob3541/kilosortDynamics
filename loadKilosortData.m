%% loadKilosortdata

% home
myKsDir = '/home/pierre/Code/Data/Neuropixel/VINNIE_COLGRID_DLPFC_NPIX45_08112023_g0/VINNIE_COLGRID_DLPFC_NPIX45_08112023_g0_imec0';
cd(myKsDir)
sp = loadKSdir(myKsDir); % from spikes repository


%work
myKsDir = '/media/pierreb/Biggie/NeuroPixel/';

% LASER/plexon data
longDirName= '/media/pierreb/Biggie/2023-09-18_09-31-11/Record Node 101/experiment1/recording1/continuous/Acquisition_Board-100.Rhythm Data/kilosort3'; 
myKsDir=longDirName;
cd(myKsDir)
rootZ = '/media/pierreb/Biggie/2023-09-18_09-31-11/Record Node 101/experiment1/recording1/continuous/Acquisition_Board-100.Rhythm Data'; % the raw data binary file is in this folder'

cd(rootZ)
sp = loadKSdir(myKsDir); % from spikes repository




%% calculate quality metrics - to determine which units waveforms to plot

% homeCodeDir
cd '/home/pierre/Code/kilosortDynamics'


% workCodeDir 
cd '~/Documents/Code/kilosortDynamics'

[spTimes,clusIdx,spTemp,sr,allChs,qualMet]=calcQuality(sp);

%% inspect waveforms of units above a threshhold of spikes
% plot units with most spikes
sortSpClust=sortrows(qualMet.nSpClus,2,'descend');

% no 0-spike templates
sortSpClust=sortSpClust(sortSpClust(:,2) ~= 0,:);

% templates with more spikes than are being extracted
nSpikes=10000;
sortSpClust=sortSpClust(sortSpClust(:,2) > nSpikes);

% templates that have some minimum FR ...

% improve algorithm to only get waveforms from best channels so you can
% extract more and quicker
for i =1:size(sortSpClust,1)
clust=sortSpClust(i);
%[gwfparams,wF(i)]=inspectWaveforms(clust,myKsDir,sp,sr,nSpikes,qualMet);
% [gwfparams,wF(i)]=inspectWaveforms(clust,myKsDir,sp,sr,qualMet);

% clust=sortSpClust(i);
[gwfparams,wF(i)]=inspectWaveforms(clust,myKsDir,sp,sr,nSpikes,qualMet);
end

% all extracted spikes with their mean from best peak2peak channel
for i =1:size(sortSpClust,1)
figure;
hold on;
plot(wF(i).bestCh');
plot(mean(wF(i).bestCh),'k','LineWidth',3);
title(['Clust: ', num2str(sortSpClust(i)), ', Ch: ', num2str(wF(i).peak2peak(1)), ', nSps: ', num2str(nSpikes)])
end

% plot all ducks in a row
plot(wF.bestCh(1,:))

wfWin=gwfparams.wfWin;

% do artifact removal before this
nKeepSpikes=100;
train=zeros(1,wf.spikeTimeKeeps(nKeepSpikes)+wfWin(2));

% train=zeros(wf.spikeTimeKeeps(nKeepSpikes)+wfWin(2),numel(wF.chs));

for i = 1: nKeepSpikes%length(spikeTimes)
    bfSp=wf.spikeTimeKeeps(i)+wfWin(1);
    afSp=wf.spikeTimeKeeps(i)+wfWin(2);

    train(1,bfSp:afSp)=wF.bestCh(i,:);
%     train(bfSp:afSp,:)=squeeze(wF.waveForms(i,:,:))';
end
offset=(numel(wF.chs):-1:1)'*100;
plot(train+offset')
plot(train)
% when are the spike times for this cluster?


% average of 10,100, and 1000 spikes
clear avgSpike
n=10;
rngs=1:n:nSpikes+1;

for i =1:length(rngs)-1
avgSpike(:,i)=mean(wF.bestCh(rngs(i):rngs(i+1)-1,:))';
end
plot(avgSpike);

clear avgResid
nSpikes=qualMet.nSpClus(clust+1,2);
nSpToAvg=5;
rngs=1:nSpToAvg:nSpikes+1;

for i =1:length(rngs)-1
avgResid(:,i)=mean(residuals(rngs(i):rngs(i+1)-1,:))';
end
plot(avgResid);

outlier_logic=isoutlier(mean(residuals,2));
plot(wF.bestCh(outlier_logic,:));
std(residuals,0,2)

% plot cluster template
exTemp=squeeze(spTemp(clust+1,:,:));

% will plot the template from non-zero channels
figure;
hold on;
%plot(exTemp(:,wF(1).peak2peak(1:10,2)));
%figure;
for i =1:20
    subplot(5,4,i);plot(exTemp(:,qualMet.filtChs(i,clust+1)));
end
plot(exTemp)

% plot cluster function
plotClustTemp(sp,qualMet,trialLogic,trialSecs)


%%
% greater than 0
qualMet.nSpClus(:,2)>0
totTime=rez.ops.sampsToRead/(sr);

% spikes/s over the experiment
expFR=qualMet.nSpClus(:,2)./totTime; 
sum(expFR>.5)


% spikes over the time that stimuli are presented
for s = 1:size(TFeventStruct,2) 
events.trialStart(s,1)=TFeventStruct(s).TrialEventTimes.TrialStart;
events.trialEnd(s,1)=TFeventStruct(s).TrialEventTimes.TrialEnd;
events.checkOn(s,1)=TFeventStruct(s).TrialEventTimes.CheckerboardDrawnTime;
events.targetOn(s,1)=TFeventStruct(s).TrialEventTimes.TargetsDrawnTime;
events.handMove(s,1)=TFeventStruct(s).TrialEventTimes.TimeOfHandMove;
end

% visualize where trials are in the recording
tStart = events.trialStart/sr;
tEnd = events.trialEnd/sr;
line([tStart tStart], [0 1])
hold on;
line([totTime totTime], [0 1],'LineStyle','--')
line([0 0], [0 1],'LineStyle','--')
axis([-100 totTime+10 -1 2])

% basically use this range for spike time logic
trialSecs=tEnd(end) - tStart(1);
trialSecs/60 % minutes with trials

trialLogic=spTimeClus > tStart(1) & spTimeClus < tEnd(end);
frLogic=(sum(trialLogic)/trialSecs)>10;

% distribution of spikes throughout the experiment
for i =1:10
 figure;
tmpSpTimeClus=spTimeClus(:,i);
tmpLogic=~isnan(spTimeClus(:,i));
y=ones(sum(tmpLogic),i);
% total amount of time in seconds
scatter(tmpSpTimeClus(tmpLogic),y)
histogram(tmpSpTimeClus(tmpLogic))
end

plotClustTemp(sp,qualMet,trialLogic,trialSecs)




%%



%plot waveform by channel
timePts=0:size(wF.waveForms,3)-1;
for chan= 1:size(chs,1)
figure;
plot(timePts,chWfMean(chan,:));
title(['channel: ', num2str(chs(chan))])
end


%plot waveforms by spikes
chChans=peak2peak(1:4,1);
for sps= 1:gwfparams.nWf %nSpikes 
figure;
plot(timePts,squeeze(waveForms(sps,chChans,:)));
end

nSpikes=5;
figure;
for i=1:size(chChans,1)
subplot(2,2,i); plot(squeeze(waveForms(1:nSpikes,chChans(i),:))');
%subplot(2,2,i); plot(mean(squeeze(waveForms(1:nSpikes,chChans(i),:))));
title(['channel: ', num2str(chChans(i))])
end






%%
% plot cluster template
exTemp=squeeze(spTemp(clust+1,:,:));

% will plot the template from non-zero channels
figure;
hold on;
plot(exTemp(:,chs));

% the cluster of spikes seems to move significantly over time
figure; 
imagesc(squeeze(wf.waveFormsMean))
set(gca, 'YDir', 'normal'); xlabel('time (samples)'); ylabel('channel number'); 
colormap(colormap_BlueWhiteRed); caxis([-1 1]*max(abs(caxis()))/2); box off;



%% get spike PSTHs
% align spikes around trial start times
% basically you have 600 trials and start; stop and checkboard times
% you want spikes around the checkerboard
% -300 ms to 700 ms after
% I have code for this that I can repurpose and use, this is where
% kilosortDynamics really starts to kick the fuck in! 

% need the spikes between 
checkOn=events.checkOn./sr;
checkOn_pre=checkOn-.3;
checkOn_post=checkOn+.7;

periCheck=[checkOn_pre checkOn_post];



nTrials=size(events.trialStart,1); 
% 2nd unit
% get spikes for all trials
unit2=spTimeClus(:,2);
for tr = 1:nTrials
indTrialLogic=unit2 > periCheck(tr,1) & unit2 < periCheck(tr,2);
spikeTimes = unit2(indTrialLogic);
spkIdx=round((spikeTimes-periCheck(tr,1))*1000);
emp=zeros(1000,1);
emp(spkIdx)=1;
spikeTimesPeriCheck(:,tr)=emp;
end
%sum(indTrialLogic) 


% convolution kernal
kb=0.03; 
mult=3;
kernal = normpdf(-mult*kb:0.001:mult*kb,0,kb); % 2* or 3* kb will capture the full normal distribution;

for tr=1:nTrials
periCheckFR(:,tr)=conv(spikeTimesPeriCheck(:,tr),kernal,'same'); %specifies that the output is same length as input
end

% pad with NaNs 
periCheckFRNaNs=periCheckFR;
for i=1:size(periCheckFR,2)
    if (RT(i) + abs(b1(2))*1000) < ((abs(b1(2))+b2(2))*1000)
       ConvPhotoNaNs((RT(i)+abs(b1(2)*1000)):((abs(b1(2))+b2(2))*1000),i)=NaN; %specifies that the output is same length as input
    else
    end
end


% can separate trials by left and right reach first
choice=vertcat(TFeventStruct.chosenSide);
Left = choice == 1;
Right = choice == 2;

% group trials by coherence
% calculate trial coherence from cue value
U = vertcat(TFeventStruct.cue); 
R=U; % # of red squares
G=abs(R-(max(U)+min(U)));
C=(abs(R-G)./(R+G))*100; % color coherence
uniqC=unique(C);

for c = 1:length(uniqC)
tmpC=C==uniqC(c);
leftCohFR(:,c)=mean(periCheckFR(:,tmpC & Left),2);
rightCohFR(:,c)=mean(periCheckFR(:,tmpC & Right),2);
end

plot(mean(periCheckFR(:,Left),2));
hold on;
plot(mean(periCheckFR(:,Right),2));

plot(leftCohFR);
hold on;
plot(rightCohFR,'LineStyle','--');

hold on;
line([300 300], [0 7], 'LineStyle','--')









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
