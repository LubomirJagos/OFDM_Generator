Tx = comm.SDRuTransmitter(...
  'Platform','N200/N210/USRP2',...
  'IPAddress', '192.168.10.2', ...
  'CenterFrequency', 2.5e6, ...
  'InterpolationFactor', 256);
hMod = comm.DPSKModulator('BitInput',true);

for counter = 1:20
  data = randi([0 1], 30, 1);
  modSignal = step(hMod, data);
  step(Tx, modSignal);
end