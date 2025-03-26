function [Reduced_W,new_yconj] = upgrade_merge_L9(W,y_conj,y1,y2,LR)
    Reduced_W = W;
    y1_conj = y_conj(y1);
    y2_conj = y_conj(y2);
    
    Lambda_1 = LR(y1);
    Lambda_2 = LR(y2);
    
    a1 = W(1,y1);
    b1 = W(1,y1_conj);
    
    if Lambda_2 < Inf
        alpha2 = Lambda_2*(a1+b1)/(Lambda_2+1);
        beta2 = (a1+b1)/(Lambda_2 + 1);
    else
        alpha2 = a1+b1;
        beta2 = 0;
    end
    
    t = @(alpha,beta,x) (1 - x)*alpha + (x)*beta;
    for j = 1:2
        Reduced_W(j,y2) = Reduced_W(j,y2) + t(alpha2,beta2,j-1);
        Reduced_W(j,y2_conj) = Reduced_W(j,y2_conj) + t(beta2,alpha2,j-1);
    end
    
    Reduced_W(:,[y1,y_conj(y1)]) = [];

    % Updating Conjugates
    yidx = 1:length(y_conj);
    yidx(y1) = y2;
    yidx(y1_conj) = y_conj(y2);

    elimination_idx = sort([y1 y1_conj]);
    yidx(elimination_idx(1)) = 0;
    yidx(yidx > elimination_idx(1)) = yidx(yidx > elimination_idx(1)) - 1;
    yidx(elimination_idx(2)) = 0;
    yidx(yidx > elimination_idx(2)-1) = yidx(yidx > elimination_idx(2)-1) - 1;

    for i = 1:length(y_conj)
        y_conj(i) = yidx(y_conj(i));
    end

    y_conj(y_conj == 0) = [];
    new_yconj = y_conj;

end
