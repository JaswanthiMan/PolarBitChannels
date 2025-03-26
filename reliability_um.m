function [mutual_infs,reliability_seq] = reliability_um(N,p,mu)

    m = log2(N);
    n = 2^m;
    %mu = 64;
    matrixArray = zeros(2,mu,n);
    
    % Degrading Merge
    mutual_informations = zeros(1,n);
    for i = 1:n
        %Q = [0.2,0.8;0.8,0.2];
        Q = [p,1-p;1-p,p];
        Qyconj = [2,1];
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
        matrixArray(:,:,i) = Q;
        mutual_informations(i) = compute_mutualinformation(matrixArray(:,:,i));
    end
    
    [mutual_infs,reliability_seq] = sort(mutual_informations);
    %writematrix(matrixArray,'DM.txt');
end