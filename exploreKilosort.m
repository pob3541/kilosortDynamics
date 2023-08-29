
%ks.rez are the results from kilosort
results=ks.rez;

%results xc, yc are coordinates of each channel on the probe

%

spikeData=results.st0;
% similarity score
ks.rez.simScore
sort(simScore(:))

% Basically you want to get rid of the diagonal

% can I plot spikes from this?