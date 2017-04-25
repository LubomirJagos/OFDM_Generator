function sampleFreqValue = Hantek6022BE_GetSampleFreq(sampleFreq)
    if (sampleFreq >= 0 && sampleFreq <= 10) sampleFreqValue = 48e6;
    elseif (sampleFreq == 11) sampleFreqValue = 16e6;
    elseif (sampleFreq == 12) sampleFreqValue = 8e6;
    elseif (sampleFreq == 13) sampleFreqValue = 4e6;
    elseif (sampleFreq >= 14 && sampleFreq <= 24) sampleFreqValue = 1e6;
    elseif (sampleFreq == 25) sampleFreqValue = 500e3;
    elseif (sampleFreq == 26) sampleFreqValue = 200e3;
    elseif (sampleFreq == 27) sampleFreqValue = 100e3;
    else sampleFreqValue = 0;               % it shows error somewhere!
    end
end
