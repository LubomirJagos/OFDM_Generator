message = 'Live long and prosper, from the Communications System Toolbox Team at MathWorks!';
numFrames = 1e2;

% Adjustable channel parameters
EbN0dB = 12;            % Channel noise level (dB)
frequencyOffset = 1e4;  % Frequency offset (Hz)
phaseOffset = 15;       % Phase offset (Degrees)
delay = 80;             % Initial sample offset for entire data stream (samples)

% Display recovered messages
displayRecoveredMsg = false;

% Enable scope visualizations
useScopes = true;

% Check for MATLAB Coder license
useCodegen = checkCodegenLicense;
if useCodegen
  fprintf(['--MATLAB Coder license found. ',...
    'Transmitter and receiver functions will be compiled for ',...
    'additional simulation acceleration.--\n']);
end

% By default the transmitter and receiver functions will be recompiled
% between every run, which is not always necessary.  To disable receiver
% compilation, change "compileIt" to false.
compileIt = useCodegen;
