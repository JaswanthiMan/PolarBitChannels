function [upgrade_merge_fin,uyconj_fin] = upgradeMergefin(Q,yconj,mu)
    matrixSize = size(Q);
    while(matrixSize(2) > mu)
        [Q,yconj] = upgradeMerge(Q,yconj);
        matrixSize = size(Q);
    end
    upgrade_merge_fin = Q;
    uyconj_fin = yconj;
end

