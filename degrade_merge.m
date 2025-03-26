function [Reduced_matrix,new_yconj] = degrade_merge(W,y_conj)
    % Degrading merge (Algorithm C)
    % Reduces the matrix size by 2
    %size_L = size(W);
    %doub_L = size_L(2);
    
    
    %finding the likelihoods
    LR = W(1,:)./W(2,:);
    
    [sorted_LR, sorted_idx] = sort(LR);
    
    %finding the L indices with LR >= 1 such that the ...
    % conjugates are not involved in the list
    y_L_idx = find_L_indices(sorted_LR,sorted_idx,y_conj);
    
    L = length(y_L_idx);
    
    cfun = @(a,b) -(a+b)*log2((a+b)/2) + a*log2(a) + b*log2(b);
    
    capacities = zeros(1,L-1);
    for i = 1:L-1
        a = W(1,y_L_idx(i)); b = W(1,y_conj(y_L_idx(i)));
        a1 = W(1,y_L_idx(i+1)); b1 = W(1,y_conj(y_L_idx(i+1)));
        aplus = a + a1; bplus = b +b1;
        capacities(i) = cfun(a,b) + cfun(a1,b1) - cfun(aplus,bplus);
    end
    
    [~, max_cap_idx] = max(capacities);
    
    yi = y_L_idx(max_cap_idx); %yi that led to the maximum capacity
    yi_conj = y_conj(yi);
    yiplus1 = y_L_idx(max_cap_idx+1);
    yiconjplus1 = y_conj(yiplus1);
    
    % Reduce the matrix size by 2
    Reduced_matrix = W;

    Reduced_matrix(:,yi) = Reduced_matrix(:,yi) + Reduced_matrix(:,yiplus1);
    Reduced_matrix(:,yi_conj) = Reduced_matrix(:,yi_conj) + Reduced_matrix(:,yiconjplus1);
    Reduced_matrix(:,[yiplus1,yiconjplus1]) = [];
    
    % Updating Conjugates
    yidx = 1:length(y_conj);
    yidx(yiplus1) = yi;
    yidx(yiconjplus1) = y_conj(yi);

    elimination_idx = sort([yiplus1 yiconjplus1]);
    yidx(elimination_idx(1)) = 0;
    yidx(yidx > elimination_idx(1)) = yidx(yidx > elimination_idx(1)) - 1;
    yidx(elimination_idx(2)) = 0;
    yidx(yidx > elimination_idx(2)-1) = yidx(yidx > elimination_idx(2)-1) - 1;

    for i = 1:length(y_conj)
        y_conj(i) = yidx(y_conj(i));
    end

    y_conj(y_conj == 0) = [];

    %disp(y_conj);

    new_yconj = y_conj;
end
