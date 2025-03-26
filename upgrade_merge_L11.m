function [Reduced_W,new_yconj] = upgrade_merge_L11(W,y_conj,LR,y_L_idx)
    %finding the likelihoods
    %LR = W(1,:)./W(2,:);
    
    %[sorted_LR, sorted_idx] = sort(LR);
    
    %finding the L indices with LR >= 1 such that the ...
    % conjugates are not involved in the list
    %y_L_idx = find_L_indices(sorted_LR,sorted_idx,y_conj);
    
    %epsilon = 1e-3;
    
    L = length(y_L_idx);

    
    diff_capacities = zeros(1,L-2);
    
    for i = 1:L-2
        Lambda_1 = LR(y_L_idx(i));
        Lambda_2 = LR(y_L_idx(i+1));
        Lambda_3 = LR(y_L_idx(i+2));
    
        pi_2 = W(1,y_L_idx(i+1)) + W(2,y_L_idx(i+1));
    
        if Lambda_3 < Inf
            term1 = (Lambda_3 - Lambda_2)*(Lambda_1 * log2(1 + (1/Lambda_1)) + log2(1 + Lambda_1));
            term2 = (Lambda_2 - Lambda_1)*(Lambda_3 * log2(1 + (1/Lambda_3)) + log2(1 + Lambda_3));
            term3 = (Lambda_1 - Lambda_3)*(Lambda_2 * log2(1 + (1/Lambda_2)) + log2(1 + Lambda_2));
            
            diff_capacities(i) = (pi_2/((Lambda_2+1)*(Lambda_1-Lambda_3))) * (term1+term2+term3);
    
        else
            diff_capacities(i) = (pi_2/(Lambda_2+1)) * (-Lambda_1 * log2(1+1/Lambda_1) - log2(1 + Lambda_1) + Lambda_2 * log2(1+1/Lambda_2) + log2(1+Lambda_2));
        end
    end
    
    [~, min_cap_idx] = min(diff_capacities);
    
    y1 = y_L_idx(min_cap_idx);
    y2 = y_L_idx(min_cap_idx +1);
    y3 = y_L_idx(min_cap_idx +2);

    
    Lambda_1 = LR(y1);
    Lambda_2 = LR(y2);
    Lambda_3 = LR(y3);

    a2 = W(1,y2);
    b2 = W(1,y_conj(y2));
    
    if Lambda_3 < Inf
        alpha_1 = Lambda_1 * (Lambda_3 * b2 - a2)/(Lambda_3 - Lambda_1);
        beta_1 = (Lambda_3 * b2 - a2)/(Lambda_3 - Lambda_1);
        alpha_3 = Lambda_3 * (a2 - Lambda_1 * b2)/(Lambda_3 - Lambda_1);
        beta_3 = (a2 - Lambda_1 * b2)/(Lambda_3 - Lambda_1);
    else
        alpha_1 = Lambda_1 * b2;
        beta_1 = b2;
        alpha_3 = a2 - Lambda_1 * b2;
        beta_3 = 0;
    end
    
    t = @(alpha,beta,x) (1 - x)*alpha + (x)*beta;
    
    Reduced_W = W;
    
    for i = 1:2
        Reduced_W(i,y1) = Reduced_W(i,y1) + t(alpha_1,beta_1,i-1);
        Reduced_W(i,y_conj(y1)) = Reduced_W(i,y_conj(y1)) + t(beta_1,alpha_1,i-1);
        Reduced_W(i,y3) = Reduced_W(i,y3) + t(alpha_3,beta_3,i-1);
        Reduced_W(i,y_conj(y3)) = Reduced_W(i,y_conj(y3)) + t(beta_3,alpha_3,i-1);
    end
    
    Reduced_W(:,[y2,y_conj(y2)]) = [];

    % Updating Conjugates
    yidx = 1:length(y_conj);
    yidx(y2) = y1;
    y2_conj = y_conj(y2);
    yidx(y2_conj) = y_conj(y1);

    elimination_idx = sort([y2 y2_conj]);
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

