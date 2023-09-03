function plotClustTemp(sp,qualMet, trialLogic,trialSecs)%,varargin)

% which channels are these clusters on 
% will want to limit plotting to good template clusters
% quality metrics 
%   1. nSpClus > 0
% nClusters=qualMet.nClusts;
% firstClus=1;
% assignopts(who,varargin)

% plot templates that have 0 spikes
plotLogic=qualMet.nSpClus(:,2) == 0;
toPlot=find(plotLogic); 

% plot templates above some estimated FR
frThresh=5; % n spikes/s from start of first trial to end of last one
frLogic=(sum(trialLogic)/trialSecs)>frThresh;
toPlot=find(frLogic); 


for clus = 1 : length(toPlot)
exTemp=squeeze(sp.temps(toPlot(clus),:,:));

% finds the nonzero channels to plot
[rows,cols] = find(exTemp);
filtChs=unique(cols);

% will plot the template from non-zero channels
figure;
hold on;
plot(exTemp(:,filtChs));
end

end