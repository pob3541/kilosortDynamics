function plotClustTemp(sp,qualMet)%,varargin)

% which channels are these clusters on 
% will want to limit plotting to good template clusters
% quality metrics 
%   1. nSpClus > 0
% nClusters=qualMet.nClusts;
% firstClus=1;
% assignopts(who,varargin)

% use logic to choose which clusters to plot
% for fun lets only plot units with 0 spikes
% can filter out templates with 0 spikes! 
% 10% of the first data set :)
sum(qualMet.nSpClus(:,2) == 0)
plotLogic=qualMet.nSpClus(:,2) == 0;
toPlot=find(plotLogic); % The first template it plotted had biological spikes though ...
% don't necessarily want to not plot it or accidently ignore it. 
% So can set a spike threshold based on amount of time in the experiment

for clus = 1 : length(toPlot)
exTemp=squeeze(sp.temps(toPlot(clus),:,:));

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

% will plot the template from non-zero channels
figure;
hold on;
for i =1:length(filtChs)
plot(sp.temps(clus,:,filtChs(i)))
end
end
end