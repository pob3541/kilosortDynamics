%% loadKilosortdata

% home
myKsDir = '/home/pierre/Code/Data/Neuropixel/VINNIE_COLGRID_DLPFC_NPIX45_08112023_g0/VINNIE_COLGRID_DLPFC_NPIX45_08112023_g0_imec0';
cd(myKsDir)
sp = loadKSdir(myKsDir); % from spikes repository


%work
myKsDir = '/media/pierreb/Biggie/NeuroPixel/';
cd(myKsDir)
sp = loadKSdir(myKsDir); % from spikes repository




%% calculate quality metrics - to determine which units waveforms to plot

% homeCodeDir
cd '/home/pierre/Code/kilosortDynamics'


% workCodeDir 
cd '~/Documents/Code/kilosortDynamics'

[spTimes,clusIdx,spTemp,sr,allChs,qualMet]=calcQuality(sp);


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
% load in waveforms 
% added floor to nSamp in getWaveForms
% which cluster to get?
clust=1;

gwfparams.dataDir = myKsDir;    % KiloSort/Phy output folder
apD = dir(fullfile(myKsDir, '*ap*.bin')); % AP band file from spikeGLX specifically
gwfparams.fileName = apD(1).name;         % .dat file containing the raw 
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 385;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
gwfparams.spikeTimes = ceil(sp.st(sp.clu==clust)*sr); % Vector of cluster spike times (in samples) same length as .spikeClusters
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
chs = allChs(~isnan(qualMet.filtChs(:,clust)));
nSpikes=size(gwfparams.spikeTimes,1);
figure;
    figure;
for sp= 113 %nSpikes 
figure;
plot(0:81,squeeze(waveForms(sp,chs([14 16]),:)));
axis([0 90 -70 40])
end

% plot cluster template
exTemp=squeeze(spTemp(clust+1,:,:));

% will plot the template from non-zero channels
figure;
hold on;
plot(exTemp(:,chs));

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
