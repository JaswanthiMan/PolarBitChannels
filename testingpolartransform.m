clc
clear all;
close all;

m = log2(8);
N = 2^m;
mu = 64;
matrixSizes = zeros(N,2,log2(N)+1);
    
for i = 1:log2(N)+1
    for j = 1:N
        matrixSizes(j, 1, i) = 2;               
        matrixSizes(j, 2, i) = 2^(2^(i-1)) * (2 ^ mod(j-1, 2^(i-1)));
    end
end

matrixArray = cell(N, log2(N)+1);

for j = 1:log2(N)+1
    for i = 1:N
    matrixArray{i,j} = zeros(matrixSizes(i,1,j), matrixSizes(i,2,j));
    end
end

%% Degrading Merge
for i = 1:N
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
        Q = W;
        Qyconj = Wyconj;
        %disp([size(Q)]);
    end
    %disp(i);
    matrixArray{i,log2(N)+1} = W;
end

matrixArray1 = polar_transformation([0.2,0.8;0.8,0.2],8);

%writematrix(matrixArray,'DM.txt');

%matrixread = reshape(readmatrix('DM.txt'),2,mu,n);