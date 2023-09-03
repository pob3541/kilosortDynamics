function [spTimes,clusIdx,spTemp,qualMet]=calcISI(sp)

spTimes=sp.st;
clusIdx = sp.clu; % same as sp.spikeTemplates I think
spTemp = sp.temps; 
qualMet.nClusts = size(spTemp,1); % will not get clusters that have 0 spikes associated

% sum number of instances of each cluster
for i =0 : max(clusIdx)
clusLogic(:,i+1)=clusIdx==i;
end
nSpClus=sum(clusLogic)';


% spike times organized by cluster
nNans=max(sum(clusLogic));
for i = 1: size(clusLogic,2)
tmpSpClus=spTimes(clusLogic(:,i));
nanFill=nNans-length(tmpSpClus);
spTimeClus(:,i)=[tmpSpClus;nan(nanFill,1)];
end
qualMet.spTimeClus=spTimeClus;


%% ISI calculations
% the rule of thumb is that if of the ISIs less than 1.5% violate the ISI
% (less  than 1.5 ms) then it is labelled as a single unit

% %sr= sp.sample_rate;

for i = 1:size(spTimeClus,2)
tmpSpClus=spTimeClus(:,i);
tmpISI=diff(tmpSpClus(~isnan(spTimeClus(:,i))))*1000; % in ms
nanFill=nNans-length(tmpISI);
ISI(:,i)=[tmpISI;nan(nanFill,1)];
end



% I want unique per column
for i = 1:size(spTimeClus,2)
unitISIs=ISI(:,i);
unitISIs=unitISIs(~isnan(unitISIs));
uniqueISIs=unique(unitISIs);
pUniqueISIs(i,1)=(length(uniqueISIs)/length(unitISIs))*100;
end

% or conversely perhaps a better metric is num times an ISI repeats


%% make this into a table Cols: 1. unit #, 2. nSpikes/Cluster, 
% 3. ISI [min median average max pISI_violation nUniqueISI/nISIs], 
% 4. spike duration / percentage of recording clust is present
% 5. average/median nSpikes/trial
qualMet.nSpClus = [(0:length(nSpClus)-1)',nSpClus];
pISIviol = (sum(ISI<1.5)'./nSpClus)*100;
varNames={'minISI','medianISI','maxISI','pISIviol'};
arrISI=[min(ISI)',median(ISI,'omitnan')',max(ISI)',pISIviol];
tblISI=array2table(arrISI,"VariableNames",varNames);

qualMet.tblISI= tblISI;



%% for plotting ISI
% figure;
%boxplot(ISI);
% histogram(ISI(ISI<1000),'Normalization','probability')
% boxplot(ISI(ISI<10))


%histogram(ISI); xlabel('Time (ms)');
%ylabel ('ISI count')

% nan logic to store all the ISIs - could do a boxplot
end
% tmpStruct=[];
