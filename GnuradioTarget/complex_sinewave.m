% Citanie komplexneho signalu zo suboru z gnuradia.

f = fopen('_complex_sinewave.txt','r');
sig = fread(f,128,'float32');
fclose(f);

% rozdelenie cisla na realnu a imaginarnu cast
sig = reshape(sig,2,length(sig)/2);

figure;
plot(sig(1,:));
hold on;
plot(sig(2,:),'-r');
hold off;






