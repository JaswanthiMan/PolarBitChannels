clc
clear all;
close all;

W = [0.2,0.8;0.8,0.2];
yconj = [2,1];

[W1,yconj1] = squareStar1(W,yconj);
[DW,Dyconj] = degrade_merge(W1,yconj1);
