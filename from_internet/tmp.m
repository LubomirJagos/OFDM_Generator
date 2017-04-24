%LuboJ.
%there is no synchronization, so I could be on any sample for beginning
%no iteration over whole samples there is same samples left in the end,
%to be sure not exceed matrix dimension

for i=1:samplesPerBit:length(demodulated)-samplesPerBit
 y=[y trapz(tc(i:i+samplesPerBit-1),demodulated(i:i+samplesPerBit-1))];
end
received=y>0;
figure;
stem(received)
title('Impulses of Received bits');
xlabel('Time (seconds)-->');
ylabel('Amplitude (volts)')

figure;
stem(b)

figure;
[sigCorr, lag] = xcorr(b,double(received));
plot(lag, sigCorr, 'k');
title('Signals correlation');

