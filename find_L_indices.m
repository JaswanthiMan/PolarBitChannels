function y_L_idx = find_L_indices(sorted_LR,sorted_idx,y_conj)
    % finding the L indiced list with out conjugates for merging
    % equation-20
    y_L_idx = [];
    l = length(sorted_LR);
    for i = 1:l
        if sorted_LR(i) >= 1
            if(~ismember(y_conj(sorted_idx(i)),y_L_idx))
                y_L_idx = [y_L_idx,sorted_idx(i)];
            end
        end
    end
end