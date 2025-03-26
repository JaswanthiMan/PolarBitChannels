function [Reduced_W,new_yconj] = upgradeMerge(W,y_conj)

    LR = W(1,:)./W(2,:);
    [sorted_LR, sorted_idx] = sort(LR);
    y_L_idx = find_L_indices(sorted_LR,sorted_idx,y_conj);
    
    L = length(y_L_idx);
    flags = zeros(1,L-1);
    for i = 1:L-1
        if (LR(y_L_idx(i+1))/LR(y_L_idx(i)) < 1+1e-3)
            flags(i) = 1;
        end
    end
    flags_idx = find(flags==1);
    
    if(sum(flags)>1)
        y1 = y_L_idx(flags_idx(1));
        y2 = y_L_idx(flags_idx(1)+1);
        [Reduced_W,new_yconj] = upgrade_merge_L9(W,y_conj,y1,y2,LR);
        %disp('lemma9');
    else
        [Reduced_W,new_yconj] = upgrade_merge_L11(W,y_conj,LR,y_L_idx);
        %disp('lemma11');
    end
end