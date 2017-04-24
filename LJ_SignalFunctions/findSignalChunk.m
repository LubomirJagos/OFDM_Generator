%% These are auxiliary functions for signal processing
%
%   Better have them on one place :)
%

function [startIndex sigCorr lag] = findSignalChunk(originalSignal, signalChunk)
    % Computing signals correlation
    [sigCorr, lag] = xcorr(originalSignal, signalChunk);
    % Magic computation of start index
    startIndex = ceil(length(lag)/2 - find(max(sigCorr) == sigCorr));
end

