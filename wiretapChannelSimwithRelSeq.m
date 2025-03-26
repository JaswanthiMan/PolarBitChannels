clc
clear all;
close all;

N = 256;
%p1 = 0.05;
p2 = 0.4;
mu = 64;
% I = mutual Information, Pe = errorProb, Z = Bhattacharyya Parameter
%[IBob,PeBob,ZBob] =computeEZI_dm(N,p1,mu);
[IEve,PeEve,ZEve] =computeEZI_um(N,p2,mu);

%writematrix([IBob;PeBob;ZBob],'IPZ_Bob_005_256.txt');
writematrix([IEve;PeEve;ZEve],'IPZ_Eve_04_256.txt');




