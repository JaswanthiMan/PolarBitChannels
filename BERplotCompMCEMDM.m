clc
clear all;
close all;

% N = 256
p_range = logspace(log10(0.05),log10(0.1),7);

BER_em = [0.0135,0.0211,0.0354,0.0680,0.1007,0.1432,0.2049];
FER_em = [0.0950,0.1220,0.2060,0.3430,0.4720,0.5920,0.7390];

BER_mc = [0.0094,0.0119,0.0253,0.0536,0.0963,0.1401,0.2075];
FER_mc = [0.0490,0.0610,0.1220,0.2140,0.3600,0.5330,0.7210];

BER_dm = [0.03,0.0397,0.0645,0.0833,0.1611,0.1830,0.2662];
FER_dm = [0.1830,0.2570,0.3670,0.4350,0.6490,0.7210,0.7580];

figure;
loglog(p_range,FER_em,'*-',LineWidth=2);
hold on;
loglog(p_range,FER_mc,'*-',LineWidth=2);
hold on;
loglog(p_range,FER_dm,'*-',LineWidth=2);
hold off;
xlabel('p');
ylabel('FER')
title('FER at R = 0.5, N = 256');
grid on;
legend('Erasure Method', 'Monte-Carlo', 'Merging');