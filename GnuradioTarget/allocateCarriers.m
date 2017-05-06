function frame = allocateCarriers(  ...
    data,                         ...
    fftLen,                         ...
    packetLen,                      ...
    nProcessPackets,                ...
    occupiedCarriers,               ...
    pilotCarriers,                  ...
    pilotSymbols,                   ...
    sync1,                          ...
    sync2)

% fftLen = 64;
% packetLen = 96;
% nProcessPackets = 1;
% occupiedCarriers = [39:43 45:57 59:64 2:7 9:21 23:27];  % <-------- It has to be this way, don't change order, otherwise it's not giving right results!
% pilotCarriers = [44 58 8 22];
% pilotSymbols = [10 10 10 -10];
% sync1 = [0., 0., 0., 0., 0., 0., 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 0., 0., 0., 0., 0.];
% sync2 = [0, 0, 0, 0, 0, 0, -1, -1, -1, -1, 1, 1, -1, -1, -1, 1, -1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, 0, 1, -1, 1, 1, 1, -1, 1, 1, 1, -1, 1, 1, 1, 1, -1, 1, -1, -1, -1, 1, -1, 1, -1, -1, -1, -1, 0, 0, 0, 0, 0];

nCarriers = length(occupiedCarriers);
nSymbols = ceil(packetLen/nCarriers);

readIndex = 1;
frame = [sync1 sync2];
dataEnd = false;
while (nProcessPackets > 0 && ~dataEnd)
    needToRead = packetLen;
    while (needToRead > 0 && ~dataEnd)
        
        stopIndex = 0;
        if (needToRead > nCarriers)
            stopIndex = readIndex + nCarriers - 1;
            count = nCarriers;
        else
            stopIndex = readIndex + needToRead - 1;
            count = needToRead;
        end
        if (stopIndex <= length(data))
            dataIn = data(readIndex:stopIndex);
            readIndex = stopIndex + 1;
        else
            dataIn = data(readIndex:end);
            count = length(data) - readIndex;
            
            dataEnd = true;
            needToRead = 0;
        end

        if (count == nCarriers)
            needToRead = needToRead - count;
        else
            dataIn = [dataIn zeros(1,nCarriers-count)];
            needToRead = 0;
        end
        fftFrame = zeros(1,fftLen);
        fftFrame(occupiedCarriers) = dataIn;
        fftFrame(pilotCarriers) = pilotSymbols;
        fftFrame = circshift(fftFrame', floor(fftLen/2))';
        frame = [frame fftFrame];
    end
    frame = [frame sync1 sync2];
    nProcessPackets = nProcessPackets-1;
end

end
