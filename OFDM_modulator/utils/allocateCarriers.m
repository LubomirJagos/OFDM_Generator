function frame = allocateCarriers(  ...
    data,                           ...
    fftLen,                         ...
    packetLen,                      ...
    nProcessPackets,                ...
    occupiedCarriers,               ...
    pilotCarriers,                  ...
    pilotSymbols,                   ...
    sync1,                          ...
    sync2)

nCarriers = length(occupiedCarriers);

frame = [];
readIndex = 1;
dataEnd = false;
while (nProcessPackets > 0 && ~dataEnd)
    frame = [frame sync1 sync2];
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
        if (~isempty(pilotCarriers))
            fftFrame(pilotCarriers) = pilotSymbols;
        end
        fftFrame = circshift(fftFrame', floor(fftLen/2))';
        frame = [frame fftFrame];
    end
    nProcessPackets = nProcessPackets-1;
end

end
