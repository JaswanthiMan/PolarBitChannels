clc
clear all;
close all;

m = log2(256);
n = 2^m;
mu = 64;
matrixArray = zeros(2,mu,n);

%% Degrading Merge
for i = 1:n
    Q = [0.2,0.8;0.8,0.2];
    Qyconj = [2,1];
    %[Q,Qyconj] = degrade_merge_fin(Q,Qyconj,mu);
    ibin = dec2bin(i-1,m);
    for j = 1:length(ibin)
        if ibin(j) == '0'
            [W,Wyconj] = squareStar1(Q,Qyconj);
        else
            [W,Wyconj] = circleStar1(Q,Qyconj);
        end
        [Q,Qyconj] = degrade_merge_fin(W,Wyconj,mu);
        disp([size(Q)]);
    end
    disp(i);
    matrixArray(:,:,i) = Q;
end

writematrix(matrixArray,'DM.txt');

%matrixread = reshape(readmatrix('DM.txt'),2,mu,n);