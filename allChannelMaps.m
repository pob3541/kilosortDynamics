

%Where do the channelMaps differ?
% 384 channels
phase3A.chanMap=chanMap;
phase3A.chanMap0ind=chanMap0ind;
phase3A.connected=connected;
phase3A.shankInd=shankInd;
phase3A.xcoords=xcoords;
phase3A.ycoords=ycoords;

%276 channels
phase3A_o4.chanMap=chanMap;
phase3A_o4.chanMap0ind=chanMap0ind;
phase3A_o4.connected=connected;
phase3A_o4.shankInd=shankInd;
phase3A_o4.xcoords=xcoords;
phase3A_o4.ycoords=ycoords;

% 384 channels
phase3B1.chanMap=chanMap;
phase3B1.chanMap0ind=chanMap0ind;
phase3B1.connected=connected;
phase3B1.shankInd=shankInd;
phase3B1.xcoords=xcoords;
phase3B1.ycoords=ycoords;

% 384 channels
phase3B2.chanMap=chanMap;
phase3B2.chanMap0ind=chanMap0ind;
phase3B2.connected=connected;
phase3B2.shankInd=shankInd;
phase3B2.xcoords=xcoords;
phase3B2.ycoords=ycoords;

% same chanMap, chanMap0ind, shank index and likely same ycoords
[phase3A.chanMap,phase3B1.chanMap ,phase3B2.chanMap]
[phase3A.chanMap0ind,phase3B1.chanMap0ind ,phase3B2.chanMap0ind]
sum([phase3A.shankInd,phase3B1.shankInd ,phase3B2.shankInd])

% seem to be same ycoords
sum(diff([phase3A.ycoords,phase3B1.ycoords ,phase3B2.ycoords]))

% difference in number of connected
sum([phase3A.connected,phase3B1.connected ,phase3B2.connected])
unique([phase3A.connected,phase3B1.connected ,phase3B2.connected])


% difference in xcoords 
sum(diff([phase3A.xcoords,phase3B1.xcoords ,phase3B2.xcoords]))

sum(diff(xcoords))

%closest to phase3B1/3B2 - NP 1.0 probes are phase3B2



