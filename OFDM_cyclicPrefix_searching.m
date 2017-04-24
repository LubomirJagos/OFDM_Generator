%%   Cyclic prefix searching
%
%   Looking through random OFDM signal for cyclic prefix by using
%   correlation

clear all;
close all;

%%  Global definitions
    numCarriers = 256;
    blockSize = numCarriers;
    cpLength = ceil(0.1*blockSize);
    
    addpath('LJ_SignalFunctions');

%%  Reading file with samples
%
%   * 256 carriers
%   * 256 block size
%   * 27 cyclic prefix length


    % Data sensed by HANTEK scope
    if (~exist('Modulated'))
        load('measurements\HANTEK256carriersQPSKsampling250ns4MHzRandomSeq_25november_2.mat');
        ofdmSig = real(Modulated);
    end
    
    % Data exported from GUI are raw OFDM encoded data
%     if (~exist('ofdmSig')) ofdmSig = dlmread('OFDm256carriersQPSKsampling250ns4MHz.mat');
%     end

%     [offset sigCorr lag] = findSignalChunk(Modulated, abs(ofdmSig));    
%     figure;
%     stem(lag, sigCorr);
%     figure;
%     hold on;
%     plot(Modulated, 'b-x');
%     plot([zeros(1, offset) (max(Modulated)/max(abs(ofdmSig)))*abs(ofdmSig)], 'r-x');
%     hold off;
    
%     figure;
%     subplot(211); plot(real(ofdmSig), '-x');
%     subplot(212); plot(imag(ofdmSig), '-x');
    
    % Signal is raw OFDM signal, so by correlation can find beginning of
    % cyclic prefix and to properly synchronize with OFDM samples. This is
    % just example when I have perfectly synchronized OFDM from beginning.
%     offset = (cpLength+blockSize)*7;                  %symbol offset to not be looking for first cycluc prefix
%     sig2 = ofdmSig(offset:offset+cpLength);
%     [shiftSig2 sigCorr lag] = findSignalChunk(ofdmSig, sig2);
%     chunkPos = findSignalChunkOccurencies(ofdmSig, sig2);    
%     ofdmChunkPosSig = zeros(1, length(ofdmSig));
%     ofdmChunkPosSig(chunkPos(1):chunkPos(1)+length(sig2)-1) = sig2;
%     ofdmChunkPosSig(chunkPos(2):chunkPos(2)+length(sig2)-1) = sig2;
% 
%     figure(1);
%     subplot(211);
%     plot(lag, sigCorr, 'k');
%     title('Signals correlation, without any offset');
%     
%     figure(2);
%     hold on;
%     stem(ofdmSig, '-x');
%     stem(ofdmChunkPosSig, 'r--o');
%     title('Overlap of original signal and founded prefix')
%     hold off;

    % Here is beginning of OFDM symbol shifted. So I have to look where is
    % beginning of cyclic prefix and find it in stream.
%     signalOffset = cpLength;                                         %random offset (block size = 256)
%     offset = (cpLength+blockSize)*7;                  %some symbol offset
%     sig2 = ofdmSig(signalOffset+offset:signalOffset+offset+cpLength*5);
%     [shiftSig2 sigCorr lag] = findSignalChunk(ofdmSig, sig2);
%     sig2 = [zeros(1, abs(shiftSig2)) sig2];
%     
%     figure(1);
%     subplot(212);
%     plot(lag, sigCorr, 'k');
%     title('Signals correlation of random signal chunk');
%     
%     figure(3);
%     hold on;
%     stem(ofdmSig, '-x');
%     stem(sig2, 'r--o');
%     title('Overlap of original signal and random chunk')
%     hold off;


    % Signal is constructed in baseband and after that is modulated by IQ
    % modulator with carrier fc, HANTEK scope sees 
    cpOccur = [];
    i = 0;
    while (i < length(ofdmSig))                    %length(ofdmSig)-1        
        %skipping already found prefixes
        if (~isempty(find((cpOccur-cpLength) == i))) i = i + 2*cpLength;
        else i = i+1
        end
        
        % handbreak :D
        if (i > length(ofdmSig) - cpLength)
            break;
        end        
        
        sigChunk = ofdmSig(i:i+cpLength-1);
        sigChunkOccur = findSignalChunkOccurencies(ofdmSig, sigChunk);
        if (length(sigChunkOccur) == 2)
            i = i + cpLength;                   %skip found prefix, because also half of prefix has also huge corrlation
            cpOccur = [cpOccur sigChunkOccur];
        end
    end
    cpOccur
    
    % Plotting ofdm signal with lines showing prefixes location
    plot(abs(ofdmSig));
    hold on;
    stem(cpOccur, 0.2*ones(1,length(cpOccur)), 'r-x');





