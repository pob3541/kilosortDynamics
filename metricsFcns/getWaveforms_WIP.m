%

% get into this binary file to determine its format
fileName
fileID = fopen(fileName);
A = fread(fileID,[gwfparams.nCh nSamp],'uint16');



% memmap way
fileName = fullfile(gwfparams.dataDir,gwfparams.fileName);           
filenamestruct = dir(fileName);
dataTypeNBytes = numel(typecast(cast(0, gwfparams.dataType), 'uint8')); % determine number of bytes per sample
nSamp = filenamestruct.bytes/(gwfparams.nCh*dataTypeNBytes);  % Number of samples per channel
wfNSamples = length(gwfparams.wfWin(1):gwfparams.wfWin(end));
mmf = memmapfile(fileName, 'Format', {gwfparams.dataType, [gwfparams.nCh floor(nSamp)], 'x'})

% 3619271 structure with [nch nsamp]/structure; where's the 3 mil from?
size(mmf.Data)
mmf.Data(1).x

