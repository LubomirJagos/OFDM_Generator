%Start by creating modulator and constellation points
hQAMMod = comm.GeneralQAMModulator;     
% Setup a three point constellation
hQAMMod.Constellation = [1 0.3+i -1];
data = randi([0 2],100,1);
modData = step(hQAMMod, data);
scatterplot(modData)




