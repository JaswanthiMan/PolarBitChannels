clc
clear all;
close all;


mu = 64;
n = 256;
degraded_matrices = reshape(readmatrix('DM.txt'),2,mu,n);
upgraded_matrices = reshape(readmatrix('UM.txt'),2,mu,n);

error_probs_DM = zeros(1,n);
error_probs_UM = zeros(1,n);

for i  = 1:n
    error_probs_DM(i) = error_prob(degraded_matrices(:,:,i));
    error_probs_UM(i) = error_prob(upgraded_matrices(:,:,i));
end