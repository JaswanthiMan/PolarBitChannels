function [degrade_merge_fin,dyconj_fin] = degrade_merge_fin(Q,yconj,mu)
    matrixSize = size(Q);
    while(matrixSize(2) > mu)
        [Q,yconj] = degrade_merge(Q,yconj);
        matrixSize = size(Q);
    end
    degrade_merge_fin = Q;
    dyconj_fin = yconj;
end