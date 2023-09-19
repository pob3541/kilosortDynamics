function avgSpike(clust,nSpToAvg,qualMet,wF)

nSpikes=qualMet.nSpClus(clust+1,2);
%nSpToAvg=10;
rngs=1:nSpToAvg:nSpikes+1;

for i =1:length(rngs)-1
avgSpike(:,i)=mean(wF.bestCh(rngs(i):rngs(i+1)-1,:))';
end
plot(avgSpike);

end