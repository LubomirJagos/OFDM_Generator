%   Returns indexes of occurencies signal chunk in signal

function [index] = findSignalChunkOccurencies(originalSignal, signalChunk)
    threshold = 0.95;
    % Computing signals correlation
    [sigCorr, lag] = xcorr(originalSignal, signalChunk);
    % Magic computation of start index    
    index = find(sigCorr > max(sigCorr)*threshold);
    index = index - ceil(length(lag)/2) + 1;
end