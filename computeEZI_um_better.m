Q = [0.2,0.8;0.8,0.2];
Qyconj = [2,1];
N = 8;
matrixSizes = zeros(N,2,log2(N)+1);

for i = 1:log2(N)+1
    for j = 1:N
        matrixSizes(j, 1, i) = 2;               
        matrixSizes(j, 2, i) = min(mu,2^(2^(i-1)) * (2 ^ mod(j-1, 2^(i-1))));
    end
end

matrixArray = cell(N, log2(N)+1);
matrixArrayconj = cell(N,log2(N)+1);

for j = 1:log2(N)+1
    for i = 1:N
        matrixArray{i,j} = zeros(matrixSizes(i,1,j), matrixSizes(i,2,j));
        matrixArrayconj{i,j} = zeros(1,matrixSizes(i,2,j));
    end
end
m = log2(N);
n = 2^m;
mu = 64;


% Degrading Merge
error_probs = zeros(1,n);
mutualinformations = zeros(1,n);
battacharyaaZs = zeros(1,n);
%Q = [p,1-p;1-p,p];
matrixArray{1,1} = Q;
matrixArrayconj{1,1} = [2,1];
error_prob1 = zeros(1,N);
for columns = 2:log2(N)+1
    N1 = 2^(columns-2);
    for k = 1:N1
        %matrixSize = size(matrixArray{2*k-1,columns});
        W_N_i = matrixArray{k,columns-1};
        W_N_i_conj = matrixArrayconj{k,columns-1};
        [SS,SSconj] = squareStar1(W_N_i,W_N_i_conj);
        [matrixArray{2*k-1,columns},matrixArrayconj{2*k-1,columns}] = upgradeMergefin(SS,SSconj,mu);
        [CS,CSconj] = circleStar1(W_N_i,W_N_i_conj);
        [matrixArray{2*k,columns},matrixArrayconj{2*k,columns}] = upgradeMergefin(CS,CSconj,mu);
        if (columns == log2(N)+1)
            error_prob1(2*k-1) = error_prob(matrixArray{2*k-1,columns});
            error_prob1(2*k) = error_prob(matrixArray{2*k,columns});
        end
    end
    disp(k);
end
matrixArray1 = zeros(2,mu,n);
for i = 1:n
    %[Q,Qyconj] = degrade_merge_fin(Q,Qyconj,mu);
    ibin = dec2bin(i-1,m);
    for j = 1:length(ibin)
        if ibin(j) == '0'
            [W,Wyconj] = squareStar1(Q,Qyconj);
        else
            [W,Wyconj] = circleStar1(Q,Qyconj);
        end
        [Q,Qyconj] = upgradeMergefin(W,Wyconj,mu);
        %disp([size(Q)]);
    end
    disp(i);
    matrixArray1(:,:,i) = Q;
    %battacharyaaZs(i) = compute_Z(matrixArray1(:,:,i));
    %mutualinformations(i) = compute_mutualinformation(matrixArray1(:,:,i));
    error_probs(i) = error_prob(matrixArray1(:,:,i));
end
%[battacharyyaZ,reliability_seq] = sort(battacharyaaZs);
%writematrix(matrixArray,'DM.txt');
