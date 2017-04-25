clear all; close all;

mu = 0;
sd = 1;
ix = -5*sd:1e-3:5*sd; %covers more than 99% of the curve

iy = pdf('normal', ix, mu, sd);
iy2 = pdf('normal', ix, 0.2, sd);
iy3 = pdf('normal', ix, 0.4, sd);
iy4 = pdf('normal', ix, 0.6, sd);
iy5 = pdf('normal', ix, 0.8, sd);

iy6 = pdf('normal', ix, mu, 1);
iy7 = pdf('normal', ix, mu, 0.7);
iy8 = pdf('normal', ix, mu, 0.5);
iy9 = pdf('normal', ix, mu, 0.2);

figure;
plot(ix,iy);
title('pdf() zmena parametru mu 0.2 --> 0.2');
grid on;
hold on;
plot(ix,iy2,'-r');
plot(ix,iy3,'-g');
plot(ix,iy4,'-black');
plot(ix,iy5,'-cyan');
hold off;

figure;
plot(ix,iy6);
title('pdf() zmena parametru sd 1 --> 0.2');
grid on;
hold on;
plot(ix,iy7,'-r');
plot(ix,iy8,'-g');
plot(ix,iy9,'--black');
hold off;

